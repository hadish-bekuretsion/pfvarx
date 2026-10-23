"""Input detection and validation.

Design rule: fail closed on corrupt or mismatched inputs; warn (do not reject)
on provisional QC cutpoints. Validation results are structured so the unified
QC layer can ingest them without re-parsing logs.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

from pfvarx.errors import InputValidationError

BGZF_MAGIC = b"\x1f\x8b"
BAM_MAGIC = b"BAM\x01"
CRAM_MAGIC = b"CRAM"
FASTQ_EXT = {".fq", ".fastq"}
FASTA_EXT = {".fa", ".fasta", ".fna", ".fas"}
VCF_EXT = {".vcf", ".g.vcf", ".gvcf"}
BCF_EXT = {".bcf"}
BAM_EXT = {".bam"}
CRAM_EXT = {".cram"}


@dataclass
class ValidationIssue:
    level: str  # PASS | WARN | FAIL
    code: str
    message: str
    field: str | None = None


@dataclass
class ValidationResult:
    path: str
    role: str
    detected_type: str
    issues: list[ValidationIssue] = field(default_factory=list)
    metadata: dict = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return not any(i.level == "FAIL" for i in self.issues)

    @property
    def status(self) -> str:
        if any(i.level == "FAIL" for i in self.issues):
            return "FAIL"
        if any(i.level == "WARN" for i in self.issues):
            return "WARN"
        return "PASS"

    def add(self, level: str, code: str, message: str, field_name: str | None = None) -> None:
        self.issues.append(ValidationIssue(level, code, message, field_name))

    def as_dict(self) -> dict:
        d = asdict(self)
        d["status"] = self.status
        d["ok"] = self.ok
        return d


def _strip_compression_suffix(name: str) -> str:
    lower = name.lower()
    for suf in (".gz", ".bgz", ".bgzf"):
        if lower.endswith(suf):
            return lower[: -len(suf)]
    return lower


def _open_maybe_gzip(path: Path, binary: bool = False):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == BGZF_MAGIC:
        return gzip.open(path, "rb" if binary else "rt")
    return open(path, "rb" if binary else "rt")


def detect_type(path: str | Path) -> str:
    """Detect file type from extension and magic bytes.

    Returns one of: fastq, fasta, bam, cram, vcf, gvcf, bcf, gzip_unknown,
    unknown.
    """
    p = Path(path)
    if not p.exists():
        raise InputValidationError(f"File does not exist: {p}", cause="Missing path")
    if p.stat().st_size == 0:
        raise InputValidationError(f"File is empty: {p}", cause="Zero-byte input")
    name = _strip_compression_suffix(p.name)
    ext = Path(name).suffix.lower()

    with open(p, "rb") as fh:
        head = fh.read(4)
    if head == BAM_MAGIC:
        return "bam"
    if head == CRAM_MAGIC:
        return "cram"
    if head[:2] == BGZF_MAGIC and ext not in {".vcf", ".gvcf", ".g.vcf", ".bcf", ".fastq", ".fq", ".fasta", ".fa"}:
        # Compressed BAM is BGZF; verify by decompressing first bytes.
        try:
            with gzip.open(p, "rb") as fh:
                inner = fh.read(4)
            if inner == BAM_MAGIC:
                return "bam"
        except OSError:
            pass

    if ext in BCF_EXT:
        return "bcf"
    if ext in BAM_EXT:
        return "bam"
    if ext in CRAM_EXT:
        return "cram"
    if ext in FASTQ_EXT:
        return "fastq"
    if ext in FASTA_EXT:
        return "fasta"
    if ext in {".g.vcf", ".gvcf"} or name.lower().endswith((".g.vcf", ".gvcf")):
        return "gvcf"
    if ext in VCF_EXT:
        return "vcf"

    # Content sniffing
    try:
        with _open_maybe_gzip(p) as fh:
            first = fh.readline(4096)
            if isinstance(first, bytes):
                first = first.decode("utf-8", "replace")
            second = fh.readline(4096)
            if isinstance(second, bytes):
                second = second.decode("utf-8", "replace")
    except (OSError, UnicodeDecodeError) as exc:
        raise InputValidationError(f"Cannot read {p}: {exc}", cause="Unreadable or corrupt file") from exc

    if first.startswith("##fileformat=VCF"):
        if "##ALT=<ID=NON_REF" in first or any("NON_REF" in s for s in (first, second)):
            return "gvcf"
        return "vcf"
    if first.startswith(">"):
        return "fasta"
    if first.startswith("@") and second and not second.startswith(">"):
        return "fastq"
    return "unknown"


def validate_fasta(path: str | Path, expected_sha256: str | None = None) -> ValidationResult:
    import pysam

    p = Path(path)
    res = ValidationResult(str(p), "reference", "fasta")
    try:
        fa = pysam.FastaFile(str(p))
    except Exception as exc:
        res.add("FAIL", "FASTA_UNREADABLE", f"pysam cannot open FASTA: {exc}")
        return res
    try:
        contigs = list(fa.references)
        lengths = dict(zip(fa.references, fa.lengths, strict=True))
        if not contigs:
            res.add("FAIL", "FASTA_NO_CONTIGS", "Reference contains no sequences")
        res.metadata["contigs"] = contigs
        res.metadata["lengths"] = lengths
        res.metadata["n_contigs"] = len(contigs)
        res.metadata["total_bp"] = sum(lengths.values())
        if "3D7" not in " ".join(contigs) and not any(c.startswith("Pf3D7_") for c in contigs):
            res.add(
                "WARN",
                "REFERENCE_NOT_3D7_NAMED",
                "Contig names do not follow the Pf3D7_ convention; assembly cannot be inferred from names alone.",
            )
        # Check for duplicate contig names
        if len(set(contigs)) != len(contigs):
            res.add("FAIL", "FASTA_DUPLICATE_CONTIGS", "Duplicate contig names in FASTA")
    finally:
        fa.close()

    if expected_sha256:
        digest = _sha256(p)
        res.metadata["sha256"] = digest
        if digest != expected_sha256.lower():
            res.add(
                "FAIL",
                "REFERENCE_SHA256_MISMATCH",
                f"Reference SHA-256 {digest} does not match configured {expected_sha256}",
            )
    return res


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


_READ_NAME = re.compile(r"^[@+]?([^\s/]+)")


def _canonical_read_name(header: str) -> str:
    m = _READ_NAME.match(header.strip())
    return m.group(1) if m else header.strip()


def validate_fastq(
    path: str | Path,
    mate: str | Path | None = None,
    *,
    max_records: int = 200_000,
    min_length: int = 1,
) -> ValidationResult:
    p = Path(path)
    res = ValidationResult(str(p), "reads", "fastq")
    try:
        fh = _open_maybe_gzip(p)
    except OSError as exc:
        res.add("FAIL", "FASTQ_UNREADABLE", f"Cannot open FASTQ: {exc}")
        return res
    n = 0
    lengths: list[int] = []
    bad = 0
    first_name = None
    last_name = None
    try:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            if not qual:
                res.add("FAIL", "FASTQ_TRUNCATED", f"Truncated record near record {n + 1}")
                break
            if isinstance(header, bytes):
                header = header.decode("utf-8", "replace")
                seq = seq.decode("utf-8", "replace")
                plus = plus.decode("utf-8", "replace")
                qual = qual.decode("utf-8", "replace")
            if not header.startswith("@") or not plus.startswith("+"):
                bad += 1
            if len(seq.rstrip()) != len(qual.rstrip()):
                bad += 1
            lengths.append(len(seq.rstrip()))
            if first_name is None:
                first_name = _canonical_read_name(header)
            last_name = _canonical_read_name(header)
            n += 1
            if n >= max_records:
                break
    finally:
        fh.close()
    res.metadata.update(
        {
            "records_scanned": n,
            "malformed_records": bad,
            "min_length": min(lengths) if lengths else 0,
            "max_length": max(lengths) if lengths else 0,
            "median_length": sorted(lengths)[len(lengths) // 2] if lengths else 0,
            "first_read_name": first_name,
            "last_read_name": last_name,
        }
    )
    if n == 0:
        res.add("FAIL", "FASTQ_EMPTY", "No FASTQ records found")
    if bad:
        res.add("FAIL", "FASTQ_MALFORMED", f"{bad} malformed records detected")
    if lengths and min(lengths) < min_length:
        res.add("WARN", "FASTQ_SHORT_READS", f"Minimum read length {min(lengths)} < {min_length} bp")

    if mate is not None:
        mres = validate_fastq(mate, max_records=max_records, min_length=min_length)
        res.metadata["mate"] = mres.as_dict()
        if not mres.ok:
            res.add("FAIL", "FASTQ_MATE_INVALID", f"Mate file failed validation: {mate}")
        else:
            if (
                res.metadata.get("first_read_name")
                and mres.metadata.get("first_read_name")
                and res.metadata["first_read_name"] != mres.metadata["first_read_name"]
            ):
                res.add(
                    "FAIL",
                    "FASTQ_MATE_NAME_MISMATCH",
                    "First read names differ between mates; possible sample swap or unsynchronised files",
                )
            if n != mres.metadata.get("records_scanned"):
                res.add(
                    "WARN",
                    "FASTQ_MATE_COUNT_MISMATCH",
                    f"Mate record counts differ ({n} vs {mres.metadata.get('records_scanned')})",
                )
    return res


def validate_bam(
    path: str | Path,
    *,
    reference_fasta: str | Path | None = None,
    require_index: bool = False,
    max_reads: int = 200_000,
) -> ValidationResult:
    import pysam

    p = Path(path)
    res = ValidationResult(str(p), "alignments", "bam")
    try:
        mode = "rc" if p.suffix.lower() == ".cram" else "rb"
        af = pysam.AlignmentFile(str(p), mode, reference_filename=str(reference_fasta) if reference_fasta else None)
    except Exception as exc:
        res.add("FAIL", "BAM_UNREADABLE", f"Cannot open alignment file: {exc}")
        return res
    try:
        header = af.header.to_dict()
        res.metadata["n_references"] = len(header.get("SQ", []))
        res.metadata["references"] = [sq["SN"] for sq in header.get("SQ", [])]
        res.metadata["read_groups"] = header.get("RG", [])
        res.metadata["sort_order"] = header.get("HD", {}).get("SO", "unknown")
        res.metadata["has_index"] = af.has_index()
        if res.metadata["sort_order"] != "coordinate":
            res.add("WARN", "BAM_NOT_COORDINATE_SORTED", "BAM is not coordinate-sorted; will be sorted before use")
        if require_index and not af.has_index():
            res.add("FAIL", "BAM_NO_INDEX", "Indexed BAM required but no index found")
        elif not af.has_index():
            res.add("WARN", "BAM_NO_INDEX", "No index found; one will be created")
        if not header.get("RG"):
            res.add("WARN", "BAM_NO_READ_GROUP", "No @RG read group; sample identity cannot be confirmed from header")

        n = 0
        mapped = 0
        dup = 0
        for read in af.fetch(until_eof=True):
            n += 1
            if not read.is_unmapped:
                mapped += 1
            if read.is_duplicate:
                dup += 1
            if n >= max_reads:
                break
        res.metadata.update({"reads_scanned": n, "mapped_scanned": mapped, "duplicates_scanned": dup})
        if n == 0:
            res.add("FAIL", "BAM_EMPTY", "Alignment file contains no reads")
        elif mapped / n < 0.5:
            res.add("WARN", "BAM_LOW_MAPPED_FRACTION", f"Mapped fraction {mapped / n:.2f} < 0.50")

        if reference_fasta:
            fa = pysam.FastaFile(str(reference_fasta))
            try:
                ref_lengths = dict(zip(fa.references, fa.lengths, strict=True))
            finally:
                fa.close()
            bam_lengths = {sq["SN"]: sq["LN"] for sq in header.get("SQ", [])}
            mismatched = [
                (name, bam_lengths[name], ref_lengths.get(name))
                for name in bam_lengths
                if name in ref_lengths and bam_lengths[name] != ref_lengths[name]
            ]
            missing = [name for name in bam_lengths if name not in ref_lengths]
            res.metadata["contigs_missing_from_reference"] = missing[:20]
            res.metadata["contigs_with_length_mismatch"] = mismatched[:20]
            if mismatched:
                res.add("FAIL", "BAM_REFERENCE_LENGTH_MISMATCH", f"{len(mismatched)} contigs differ in length from reference")
            if missing:
                res.add("FAIL", "BAM_CONTIGS_NOT_IN_REFERENCE", f"{len(missing)} BAM contigs absent from reference")
    finally:
        af.close()
    return res


def validate_vcf(
    path: str | Path,
    *,
    reference_fasta: str | Path | None = None,
    check_ref: bool = True,
    max_records: int = 50_000,
) -> ValidationResult:
    import pysam

    p = Path(path)
    res = ValidationResult(str(p), "variants", "vcf")
    try:
        vf = pysam.VariantFile(str(p))
    except Exception as exc:
        res.add("FAIL", "VCF_UNREADABLE", f"Cannot open variant file: {exc}")
        return res
    try:
        header = vf.header
        res.metadata["samples"] = list(header.samples)
        res.metadata["contigs"] = list(header.contigs)
        res.metadata["fileformat"] = str(header.info.get("fileformat", "unknown")) if False else str(header)
        res.metadata["has_index"] = bool(vf.index) if getattr(vf, "index", None) else False
        if not header.samples:
            res.add("WARN", "VCF_NO_SAMPLES", "Variant file has no sample columns; site-only import")

        # VCF coordinate order follows the contig dictionary, not lexical
        # string order (for example, chr2 precedes chr10 in many dictionaries).
        contig_order = {name: i for i, name in enumerate(header.contigs)}
        fallback_order: dict[str, int] = {}

        fa = pysam.FastaFile(str(reference_fasta)) if reference_fasta and check_ref else None
        n = 0
        ref_mismatch = 0
        unsorted = 0
        prev = None
        try:
            for rec in vf:
                n += 1
                if rec.contig not in contig_order and rec.contig not in fallback_order:
                    fallback_order[rec.contig] = len(contig_order) + len(fallback_order)
                rank = contig_order.get(rec.contig, fallback_order.get(rec.contig, 0))
                key = (rank, rec.pos)
                if prev is not None and key < prev:
                    unsorted += 1
                prev = key
                if fa is not None:
                    try:
                        reference_allele = fa.fetch(
                            rec.contig, rec.pos - 1, rec.pos - 1 + len(rec.ref)
                        ).upper()
                    except (KeyError, ValueError):
                        res.add("FAIL", "VCF_CONTIG_NOT_IN_REFERENCE", f"Contig {rec.contig} absent from reference")
                        break
                    if reference_allele != rec.ref.upper():
                        ref_mismatch += 1
                if n >= max_records:
                    break
        finally:
            if fa is not None:
                fa.close()
        res.metadata.update({"records_scanned": n, "ref_mismatches": ref_mismatch, "out_of_order": unsorted})
        if n == 0:
            res.add("WARN", "VCF_EMPTY", "Variant file contains no records")
        if ref_mismatch:
            res.add("FAIL", "VCF_REF_MISMATCH", f"{ref_mismatch} records have REF alleles inconsistent with the reference")
        if unsorted:
            res.add("FAIL", "VCF_UNSORTED", f"{unsorted} records are out of coordinate order")
        if str(p).lower().endswith((".g.vcf", ".gvcf")):
            res.metadata["is_gvcf"] = True
            res.add(
                "WARN",
                "GVCF_IMPORT",
                "gVCF imported as variant-only evidence; symbolic NON_REF blocks are not alternate alleles",
            )
    finally:
        vf.close()
    return res


def validate_reference_consistency(
    reference_fasta: str | Path,
    *,
    bam: str | Path | None = None,
    vcf: str | Path | None = None,
) -> ValidationResult:
    """Cross-check contig names/lengths between reference and evidence files."""
    import pysam

    fa = pysam.FastaFile(str(reference_fasta))
    try:
        ref_lengths = dict(zip(fa.references, fa.lengths, strict=True))
    finally:
        fa.close()
    res = ValidationResult(str(reference_fasta), "reference_consistency", "fasta")
    if bam:
        b = validate_bam(bam, reference_fasta=reference_fasta, require_index=False, max_reads=1)
        for issue in b.issues:
            if issue.level == "FAIL":
                res.add("FAIL", issue.code, issue.message)
    if vcf:
        v = validate_vcf(vcf, reference_fasta=reference_fasta, check_ref=False, max_records=1)
        for issue in v.issues:
            if issue.level == "FAIL":
                res.add("FAIL", issue.code, issue.message)
    res.metadata["reference_contigs"] = len(ref_lengths)
    return res


# ---------------------------------------------------------------------------
# Samplesheet
# ---------------------------------------------------------------------------
SAMPLESHEET_COLUMNS = {
    "sample_id",
    "assay",
    "input_type",
    "fastq_1",
    "fastq_2",
    "bam",
    "vcf",
    "metadata_json",
}
REQUIRED_SAMPLESHEET_COLUMNS = {"sample_id", "input_type"}


@dataclass
class SampleRecord:
    sample_id: str
    input_type: str
    assay: str = ""
    fastq_1: str | None = None
    fastq_2: str | None = None
    bam: str | None = None
    vcf: str | None = None
    metadata: dict = field(default_factory=dict)

    def evidence_path(self) -> str | None:
        return self.fastq_1 or self.bam or self.vcf


def read_samplesheet(path: str | Path, base_dir: str | Path | None = None) -> tuple[list[SampleRecord], ValidationResult]:
    p = Path(path)
    res = ValidationResult(str(p), "samplesheet", "tsv")
    if not p.exists():
        res.add("FAIL", "SAMPLESHEET_MISSING", f"Samplesheet not found: {p}")
        return [], res
    delimiter = "\t" if p.suffix.lower() in {".tsv", ".txt"} else ","
    with p.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            res.add("FAIL", "SAMPLESHEET_EMPTY", "Samplesheet has no header")
            return [], res
        fields = {f.strip() for f in reader.fieldnames}
        missing = REQUIRED_SAMPLESHEET_COLUMNS - fields
        unknown = fields - SAMPLESHEET_COLUMNS
        if missing:
            res.add("FAIL", "SAMPLESHEET_MISSING_COLUMNS", f"Missing required columns: {sorted(missing)}")
            return [], res
        if unknown:
            res.add("WARN", "SAMPLESHEET_UNKNOWN_COLUMNS", f"Ignoring unknown columns: {sorted(unknown)}")
        base = Path(base_dir) if base_dir else p.parent
        records: list[SampleRecord] = []
        seen: set[str] = set()
        for i, row in enumerate(reader, start=2):
            sid = (row.get("sample_id") or "").strip()
            if not sid:
                res.add("FAIL", "SAMPLESHEET_MISSING_SAMPLE_ID", f"Row {i} has no sample_id")
                continue
            if sid in seen:
                res.add("FAIL", "SAMPLESHEET_DUPLICATE_SAMPLE", f"Duplicate sample_id {sid!r} at row {i}")
                continue
            seen.add(sid)
            itype = (row.get("input_type") or "auto").strip().lower()
            if itype not in {"auto", "fastq", "bam", "cram", "vcf", "bcf", "gvcf"}:
                res.add("FAIL", "SAMPLESHEET_BAD_INPUT_TYPE", f"Row {i}: invalid input_type {itype!r}")
            def _resolve(val: str | None) -> str | None:
                if not val or not val.strip():
                    return None
                q = Path(val.strip())
                return str(q if q.is_absolute() else base / q)
            meta = {}
            if row.get("metadata_json"):
                import json

                try:
                    meta = json.loads(row["metadata_json"])
                except json.JSONDecodeError:
                    res.add("WARN", "SAMPLESHEET_BAD_METADATA", f"Row {i}: metadata_json is not valid JSON")
            records.append(
                SampleRecord(
                    sample_id=sid,
                    input_type=itype,
                    assay=(row.get("assay") or "").strip(),
                    fastq_1=_resolve(row.get("fastq_1")),
                    fastq_2=_resolve(row.get("fastq_2")),
                    bam=_resolve(row.get("bam")),
                    vcf=_resolve(row.get("vcf")),
                    metadata=meta,
                )
            )
        res.metadata["n_samples"] = len(records)
        if not records:
            res.add("FAIL", "SAMPLESHEET_NO_SAMPLES", "No usable sample rows")
    return records, res
