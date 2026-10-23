"""End-to-end per-sample calling orchestration.

Ties candidate discovery, evidence extraction, statistical evaluation, FDR
control and output writing into one auditable function.  This is the code path
exercised by both the CLI and the Snakemake ``call_variants`` rule.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from pfvarx.calling.evidence import Candidate, discover_indels, discover_snps, extract_site_evidence
from pfvarx.calling.policy import SiteCall, apply_fdr, evaluate_site, summarise_calls
from pfvarx.calling.writer import write_evidence_table, write_vcf
from pfvarx.config.schema import PfVarXConfig


@dataclass
class CallingResult:
    sample_id: str
    calls: list[SiteCall] = field(default_factory=list)
    vcf_path: str | None = None
    evidence_path: str | None = None
    summary: dict = field(default_factory=dict)


def _in_callable(chrom: str, pos: int, intervals: dict[str, list[tuple[int, int]]] | None) -> bool:
    if intervals is None:
        return True
    return any(start <= pos <= end for start, end in intervals.get(chrom, []))


def _reference_matches(
    fasta_path: str | None, chrom: str, pos: int, ref: str
) -> bool:
    """Verify the candidate's REF allele against the pinned FASTA.

    Returns ``True`` when no FASTA is supplied (the caller cannot contradict
    the candidate) but records the limitation; returns ``False`` on a genuine
    mismatch so the site is failed closed with ``ReferenceMismatch``.
    """
    if not fasta_path:
        return True
    import pysam

    try:
        with pysam.FastaFile(fasta_path) as fa:
            observed = fa.fetch(chrom, pos - 1, pos - 1 + len(ref)).upper()
    except (KeyError, ValueError):
        return False
    return observed == ref.upper()


def _umi_model_valid(evidence) -> bool:
    """Whether molecule-level UMI evidence supports a validated call.

    The current implementation performs fragment-level de-duplication only and
    does not build molecule consensus; therefore UMI evidence is treated as
    *unverified* and the site fails closed rather than silently passing.
    """
    return evidence.umi_status == "UMI_CONSENSUS_VALID"


def _hypothesis_family_size(
    bam_path: str,
    cfg: PfVarXConfig,
    *,
    chrom: str | None,
    start: int | None,
    end: int | None,
) -> int:
    """Upper-bound possible SNP and configured anchored-INDEL hypotheses.

    This bound is based on the full requested reference span, before any
    evidence/candidate thresholds or ``max_candidates`` truncation. Per anchor
    it counts 3 SNP alleles plus every possible insertion/deletion sequence up
    to the configured maximum length; this intentionally overcounts.
    """
    import pysam

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        names = list(bam.references)
        lengths = dict(zip(bam.references, bam.lengths, strict=True))
    if chrom is not None and chrom not in lengths:
        raise ValueError(f"Contig {chrom!r} not present in BAM")
    scanned_bases = 0
    for name in ([chrom] if chrom else names):
        contig_length = lengths[name]
        lo = max((start or 1) - 1, 0)
        hi = min(end or contig_length, contig_length)
        scanned_bases += max(0, hi - lo)
    max_indel = int(cfg.candidate.max_small_indel_bp)
    per_anchor = 3 + 2 * sum(4**length for length in range(1, max_indel + 1))
    return max(1, scanned_bases * per_anchor)


def load_callable_mask(bed_path: str | None) -> dict[str, list[tuple[int, int]]] | None:
    if not bed_path:
        return None
    mask: dict[str, list[tuple[int, int]]] = {}
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            mask.setdefault(parts[0], []).append((int(parts[1]) + 1, int(parts[2])))
    return mask


def call_sample(
    bam_path: str,
    cfg: PfVarXConfig,
    *,
    sample_id: str,
    fasta_path: str | None = None,
    callable_bed: str | None = None,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    max_candidates: int | None = None,
    compute_bootstrap: bool | None = None,
    write_outputs: bool = True,
    output_dir: str | None = None,
) -> CallingResult:
    """Run the full calling chain for one sample."""
    if compute_bootstrap is None:
        compute_bootstrap = cfg.statistics.lrt_null_calibration == "bootstrap"

    candidates: list[Candidate] = []
    candidates.extend(
        discover_snps(
            bam_path, cfg, chrom=chrom, start=start, end=end,
            fasta_path=fasta_path, max_candidates=max_candidates,
        )
    )
    candidates.extend(
        discover_indels(
            bam_path, cfg, chrom=chrom, start=start, end=end,
            fasta_path=fasta_path, max_candidates=max_candidates,
        )
    )

    mask = load_callable_mask(callable_bed)
    calls: list[SiteCall] = []
    for cand in candidates:
        ev = extract_site_evidence(
            bam_path, cand.chrom, cand.pos, cand.ref, cand.alt, cfg,
            sample_id=sample_id, fasta_path=fasta_path,
        )
        call = evaluate_site(
            ev,
            cfg,
            in_callable_region=_in_callable(cand.chrom, cand.pos, mask),
            reference_matches=_reference_matches(fasta_path, cand.chrom, cand.pos, cand.ref),
            umi_model_valid=_umi_model_valid(ev),
            compute_bootstrap=compute_bootstrap,
        )
        calls.append(call)

    family_size = _hypothesis_family_size(
        bam_path, cfg, chrom=chrom, start=start, end=end
    )
    calls = apply_fdr(calls, cfg, family_size=family_size)
    result = CallingResult(sample_id=sample_id, calls=calls)
    result.summary = summarise_calls(calls)

    if write_outputs:
        out_dir = output_dir or cfg.outputs.root
        result.vcf_path = write_vcf(
            calls, cfg, f"{out_dir}/{sample_id}/{sample_id}.pfvarx.vcf",
            sample_id=sample_id,
        )
        result.evidence_path = write_evidence_table(
            calls, f"{out_dir}/{sample_id}/{sample_id}.evidence.tsv"
        )
    return result
