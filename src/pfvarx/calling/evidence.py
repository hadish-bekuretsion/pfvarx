"""Allele-level evidence extraction from read alignments.

Key correctness properties:

* CIGAR is parsed explicitly; soft clips are not counted as independent SNP
  evidence and insertions/deletions are evaluated at their anchor base.
* Overlapping mates from one physical fragment contribute at most one
  independent observation (``count_overlapping_mates_once``).
* Secondary/supplementary alignments are excluded from independent counts.
* Duplicate-marked reads are excluded from evidence only when the configured
  policy says so; the source BAM is never modified.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import pysam

from pfvarx.config.schema import PfVarXConfig
from pfvarx.errors import InputValidationError

BASES = "ACGT"


@dataclass
class SiteEvidence:
    sample_id: str
    chrom: str
    pos: int  # 1-based
    ref: str
    alt: str
    var_type: str
    ref_count: int = 0
    alt_count: int = 0
    other_count: int = 0
    depth: int = 0
    informative_depth: int = 0
    alt_fraction: float = 0.0
    af_ci_low: float = 0.0
    af_ci_high: float = 0.0
    mean_base_quality_alt: float = 0.0
    mean_base_quality_ref: float = 0.0
    mean_mapping_quality: float = 0.0
    alt_forward: int = 0
    alt_reverse: int = 0
    ref_forward: int = 0
    ref_reverse: int = 0
    mean_alt_read_position: float = 0.0
    mean_ref_read_position: float = 0.0
    read_position_bias: float = 0.0
    strand_bias_pvalue: float | None = None
    homopolymer_length: int = 0
    left_aligned: bool = True
    duplicate_fragments_excluded: int = 0
    low_mapq_excluded: int = 0
    low_bq_excluded: int = 0
    secondary_excluded: int = 0
    sequencing_mode: str = "wgs"
    umi_status: str = "NOT_APPLICABLE"
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)


def _passes_read_filters(record: pysam.AlignedSegment, cfg: PfVarXConfig) -> tuple[bool, str]:
    if record.is_unmapped:
        return False, "unmapped"
    if record.is_secondary or record.is_supplementary:
        return False, "secondary"
    if record.mapping_quality < cfg.preprocessing.mapping_quality_min:
        return False, "low_mapq"
    if (
        cfg.preprocessing.duplicate_exclusion_from_evidence
        and cfg.run.mode == "wgs"
        and record.is_duplicate
    ):
        return False, "duplicate"
    return True, ""


def _cigar_allele(
    record: pysam.AlignedSegment,
    pos: int,
    ref: str,
    alt: str,
) -> tuple[str, int, int] | None:
    """Classify a read as ref/alt/other for a SNP or anchored small INDEL."""
    ref = ref.upper()
    alt = alt.upper()
    var_type = "SNP" if len(ref) == len(alt) == 1 else ("INS" if len(ref) == 1 else "DEL")
    seq = record.query_sequence
    quals = record.query_qualities
    if seq is None:
        return None
    ref_pos = record.reference_start  # 0-based
    q_pos = 0
    anchor0 = pos - 1
    last_base = None
    last_bq = 0
    for op, length in record.cigartuples or []:
        if op in (0, 7, 8):  # M, =, X
            if ref_pos <= anchor0 < ref_pos + length:
                offset = anchor0 - ref_pos
                qidx = q_pos + offset
                if qidx >= len(seq):
                    return None
                base = seq[qidx].upper()
                bq = int(quals[qidx]) if quals is not None else 0
                if var_type == "SNP":
                    if base == alt:
                        return ("alt", bq, qidx)
                    if base == ref:
                        return ("ref", bq, qidx)
                    return ("other", bq, qidx)
                last_base = base
                last_bq = bq
                if base != ref[0]:
                    return ("other", bq, qidx)
            ref_pos += length
            q_pos += length
        elif op == 1:  # insertion, occurs after ref_pos-1
            if var_type == "INS" and ref_pos - 1 == anchor0 and last_base == ref[0]:
                inserted = seq[q_pos : q_pos + length].upper()
                if inserted == alt[1:]:
                    return ("alt", last_bq, q_pos - 1)
                return ("other", last_bq, q_pos - 1)
            q_pos += length
        elif op == 2:  # deletion
            if var_type == "DEL" and ref_pos == anchor0 + 1 and length == len(ref) - 1:
                if last_base == ref[0]:
                    return ("alt", last_bq, max(q_pos - 1, 0))
                return ("other", last_bq, max(q_pos - 1, 0))
            ref_pos += length
        elif op == 3:  # ref skip
            ref_pos += length
        elif op == 4:  # soft clip
            q_pos += length
        elif op in (5, 6):  # hard clip / padding
            continue
    return None


def _homopolymer_length(fasta: pysam.FastaFile | None, chrom: str, pos: int, ref: str) -> int:
    if fasta is None:
        return 0
    try:
        start = max(pos - 30, 1)
        end = pos + 30
        seq = fasta.fetch(chrom, start - 1, end).upper()
    except (KeyError, ValueError):
        return 0
    if not seq:
        return 0
    centre = pos - start
    if centre < 0 or centre >= len(seq):
        return 0
    base = seq[centre]
    left = centre
    while left > 0 and seq[left - 1] == base:
        left -= 1
    right = centre
    while right + 1 < len(seq) and seq[right + 1] == base:
        right += 1
    return right - left + 1


def _is_left_aligned(fasta: pysam.FastaFile | None, chrom: str, pos: int, ref: str, alt: str) -> bool:
    """Heuristic check that an INDEL is already left-aligned."""
    if fasta is None or len(ref) == len(alt):
        return True
    try:
        prev = fasta.fetch(chrom, pos - 2, pos - 1).upper()
    except (KeyError, ValueError):
        return True
    if not prev:
        return True
    if len(ref) > len(alt):  # deletion
        return prev != ref[-1]
    return prev != alt[-1]


def extract_site_evidence(
    bam_path: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    cfg: PfVarXConfig,
    *,
    sample_id: str = "sample",
    fasta_path: str | None = None,
) -> SiteEvidence:
    """Collect allele-level evidence at one site."""
    ref = ref.upper()
    alt = alt.upper()
    if len(ref) == len(alt) == 1:
        var_type = "SNP"
    elif len(ref) == 1 and 1 < len(alt) <= cfg.candidate.max_small_indel_bp + 1:
        var_type = "INS"
    elif len(alt) == 1 and 1 < len(ref) <= cfg.candidate.max_small_indel_bp + 1:
        var_type = "DEL"
    else:
        raise InputValidationError(
            f"Variant {chrom}:{pos} {ref}>{alt} is outside the SNP/small-INDEL scope",
            cause="Allele length outside 1..max_small_indel_bp",
        )

    ev = SiteEvidence(
        sample_id=sample_id, chrom=chrom, pos=pos, ref=ref, alt=alt, var_type=var_type,
        sequencing_mode=cfg.run.mode,
        umi_status=(
            "NOT_APPLICABLE" if cfg.run.mode != "tbs_umi" else "UMI_UNVERIFIED"
        ),
    )
    try:
        af = pysam.AlignmentFile(bam_path, "rb")
    except Exception as exc:
        raise InputValidationError(f"Cannot open BAM {bam_path}: {exc}") from exc

    fasta = pysam.FastaFile(fasta_path) if fasta_path else None
    seen_fragments: set[str] = set()
    alt_bq: list[int] = []
    ref_bq: list[int] = []
    mapqs: list[int] = []
    alt_positions: list[float] = []
    ref_positions: list[float] = []
    try:
        try:
            reads = af.fetch(chrom, max(pos - 2, 0), pos + max(len(ref), 1) + 1)
        except ValueError as exc:
            raise InputValidationError(
                f"Contig {chrom} not present in BAM {bam_path}", cause=str(exc)
            ) from exc
        for record in reads:
            ok, reason = _passes_read_filters(record, cfg)
            if not ok:
                if reason == "duplicate":
                    ev.duplicate_fragments_excluded += 1
                elif reason == "low_mapq":
                    ev.low_mapq_excluded += 1
                elif reason == "secondary":
                    ev.secondary_excluded += 1
                continue
            if cfg.preprocessing.count_overlapping_mates_once and record.query_name in seen_fragments:
                continue
            classified = _cigar_allele(record, pos, ref, alt)
            if classified is None:
                continue
            allele, bq, qidx = classified
            if bq < cfg.preprocessing.base_quality_min:
                ev.low_bq_excluded += 1
                continue
            if cfg.preprocessing.count_overlapping_mates_once:
                seen_fragments.add(record.query_name)
            mapqs.append(record.mapping_quality)
            read_len = record.query_length or 1
            rel_pos = qidx / read_len if read_len else 0.0
            if allele == "alt":
                ev.alt_count += 1
                alt_bq.append(bq)
                alt_positions.append(rel_pos)
                if record.is_reverse:
                    ev.alt_reverse += 1
                else:
                    ev.alt_forward += 1
            elif allele == "ref":
                ev.ref_count += 1
                ref_bq.append(bq)
                ref_positions.append(rel_pos)
                if record.is_reverse:
                    ev.ref_reverse += 1
                else:
                    ev.ref_forward += 1
            else:
                ev.other_count += 1
    finally:
        af.close()
        if fasta is not None:
            fasta.close()

    ev.depth = ev.ref_count + ev.alt_count + ev.other_count
    ev.informative_depth = ev.ref_count + ev.alt_count
    ev.alt_fraction = ev.alt_count / ev.informative_depth if ev.informative_depth else 0.0
    if ev.informative_depth:
        from pfvarx.models.count_models import allele_fraction_interval

        lo, hi = allele_fraction_interval(ev.alt_count, ev.informative_depth)
        ev.af_ci_low, ev.af_ci_high = lo, hi
    ev.mean_base_quality_alt = sum(alt_bq) / len(alt_bq) if alt_bq else 0.0
    ev.mean_base_quality_ref = sum(ref_bq) / len(ref_bq) if ref_bq else 0.0
    ev.mean_mapping_quality = sum(mapqs) / len(mapqs) if mapqs else 0.0
    ev.mean_alt_read_position = sum(alt_positions) / len(alt_positions) if alt_positions else 0.0
    ev.mean_ref_read_position = sum(ref_positions) / len(ref_positions) if ref_positions else 0.0
    ev.read_position_bias = abs(ev.mean_alt_read_position - ev.mean_ref_read_position)

    if ev.alt_count + ev.ref_count > 0:
        from scipy import stats as _st

        table = [[ev.alt_forward, ev.alt_reverse], [ev.ref_forward, ev.ref_reverse]]
        try:
            ev.strand_bias_pvalue = float(_st.fisher_exact(table)[1])
        except ValueError:
            ev.strand_bias_pvalue = None

    if fasta_path:
        fa = pysam.FastaFile(fasta_path)
        try:
            ev.homopolymer_length = _homopolymer_length(fa, chrom, pos, ref)
            ev.left_aligned = _is_left_aligned(fa, chrom, pos, ref, alt)
        finally:
            fa.close()
    if var_type == "SNP" and ev.alt_count == 0:
        ev.notes.append("No alternate observations")
    return ev


# ---------------------------------------------------------------------------
# Candidate discovery
# ---------------------------------------------------------------------------
@dataclass
class Candidate:
    chrom: str
    pos: int
    ref: str
    alt: str
    var_type: str
    ref_count: int
    alt_count: int
    depth: int
    source: str = "pileup"

    @property
    def alt_fraction(self) -> float:
        d = self.ref_count + self.alt_count
        return self.alt_count / d if d else 0.0

    def as_dict(self) -> dict:
        d = asdict(self)
        d["alt_fraction"] = self.alt_fraction
        return d


def discover_snps(
    bam_path: str,
    cfg: PfVarXConfig,
    *,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    fasta_path: str | None = None,
    max_candidates: int | None = None,
) -> list[Candidate]:
    """Discover SNP candidates using per-position independent fragment counts."""
    af = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path) if fasta_path else None
    out: list[Candidate] = []
    try:
        contigs = [chrom] if chrom else list(af.references)
        for ctg in contigs:
            length = af.get_reference_length(ctg)
            s = max((start or 1) - 1, 0)
            e = min(end or length, length)
            counts: dict[int, dict[str, int]] = {}
            seen: dict[int, set[str]] = {}
            for record in af.fetch(ctg, s, e):
                ok, _ = _passes_read_filters(record, cfg)
                if not ok:
                    continue
                seq = record.query_sequence
                quals = record.query_qualities
                if seq is None:
                    continue
                ref_pos = record.reference_start
                q_pos = 0
                for op, ln in record.cigartuples or []:
                    if op in (0, 7, 8):
                        for off in range(ln):
                            p0 = ref_pos + off
                            if p0 < s or p0 >= e:
                                continue
                            qidx = q_pos + off
                            if qidx >= len(seq):
                                continue
                            bq = int(quals[qidx]) if quals is not None else 0
                            if bq < cfg.preprocessing.base_quality_min:
                                continue
                            if cfg.preprocessing.count_overlapping_mates_once:
                                fs = seen.setdefault(p0, set())
                                if record.query_name in fs:
                                    continue
                                fs.add(record.query_name)
                            base = seq[qidx].upper()
                            if base not in BASES:
                                continue
                            d = counts.setdefault(p0, {})
                            d[base] = d.get(base, 0) + 1
                        ref_pos += ln
                        q_pos += ln
                    elif op == 1:
                        q_pos += ln
                    elif op == 2 or op == 3:
                        ref_pos += ln
                    elif op == 4:
                        q_pos += ln
            for p0, base_counts in counts.items():
                if fasta is not None:
                    rb = fasta.fetch(ctg, p0, p0 + 1).upper()
                else:
                    rb = max(base_counts.items(), key=lambda kv: kv[1])[0]
                if rb not in BASES:
                    continue
                ref_count = base_counts.get(rb, 0)
                for b, c in base_counts.items():
                    if b == rb:
                        continue
                    depth = ref_count + c
                    if c < cfg.candidate.snp_min_alt_observations:
                        continue
                    if depth and c / depth < cfg.candidate.exploratory_min_alt_fraction:
                        continue
                    out.append(Candidate(ctg, p0 + 1, rb, b, "SNP", ref_count, c, depth))
                    if max_candidates and len(out) >= max_candidates:
                        return out
    finally:
        af.close()
        if fasta is not None:
            fasta.close()
    return out


def discover_indels(
    bam_path: str,
    cfg: PfVarXConfig,
    *,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    fasta_path: str | None = None,
    max_candidates: int | None = None,
) -> list[Candidate]:
    """Discover anchored small-INDEL candidates from CIGAR events."""
    af = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(fasta_path) if fasta_path else None
    events: dict[tuple, set[str]] = {}
    try:
        contigs = [chrom] if chrom else list(af.references)
        for ctg in contigs:
            length = af.get_reference_length(ctg)
            s = max((start or 1) - 1, 0)
            e = min(end or length, length)
            for record in af.fetch(ctg, s, e):
                ok, _ = _passes_read_filters(record, cfg)
                if not ok:
                    continue
                seq = record.query_sequence
                if seq is None:
                    continue
                ref_pos = record.reference_start
                q_pos = 0
                for op, ln in record.cigartuples or []:
                    if op in (0, 7, 8):
                        ref_pos += ln
                        q_pos += ln
                    elif op == 1:
                        anchor0 = ref_pos - 1
                        if (
                            s <= anchor0 < e
                            and 1 <= ln <= cfg.candidate.max_small_indel_bp
                            and fasta is not None
                        ):
                            anchor = fasta.fetch(ctg, anchor0, anchor0 + 1).upper()
                            inserted = seq[q_pos : q_pos + ln].upper()
                            if anchor in BASES and set(inserted) <= set(BASES):
                                events.setdefault(
                                    (ctg, anchor0 + 1, anchor, anchor + inserted, "INS"), set()
                                ).add(record.query_name)
                        q_pos += ln
                    elif op == 2:
                        anchor0 = ref_pos - 1
                        if (
                            s <= anchor0 < e
                            and 1 <= ln <= cfg.candidate.max_small_indel_bp
                            and fasta is not None
                        ):
                            anchor = fasta.fetch(ctg, anchor0, anchor0 + 1).upper()
                            deleted = fasta.fetch(ctg, ref_pos, ref_pos + ln).upper()
                            if anchor in BASES and set(deleted) <= set(BASES):
                                events.setdefault(
                                    (ctg, anchor0 + 1, anchor + deleted, anchor, "DEL"), set()
                                ).add(record.query_name)
                        ref_pos += ln
                    elif op == 3:
                        ref_pos += ln
                    elif op == 4:
                        q_pos += ln
    finally:
        af.close()
        if fasta is not None:
            fasta.close()

    out: list[Candidate] = []
    for (ctg, pos, ref, alt, vtype), frags in events.items():
        if len(frags) < cfg.candidate.indel_min_alt_observations:
            continue
        ev = extract_site_evidence(bam_path, ctg, pos, ref, alt, cfg, fasta_path=fasta_path)
        out.append(Candidate(ctg, pos, ref, alt, vtype, ev.ref_count, ev.alt_count, ev.depth, source="cigar"))
        if max_candidates and len(out) >= max_candidates:
            break
    return out
