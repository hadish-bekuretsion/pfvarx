"""Variant-import route.

Pre-called VCF/BCF/gVCF inputs must **never** be realigned.  This module reads
such a callset with pysam and converts each record into the PfVarX evidence
schema so that it flows through the same QC, annotation and reporting layers.

Because an imported callset may not carry the raw reference/alternate
observation counts that PfVarX uses for its statistical model, the importer is
explicit about what it could and could not recover:

* If ``AD`` (allelic depths) is present, ref/alt counts are taken from it.
* Else if ``AO``/``RO`` (FreeBayes-style) are present, those are used.
* Else counts are left as ``None`` and the record is tagged
  ``IMPORTED_COUNTS_UNAVAILABLE``; such sites can be reported and annotated but
  cannot receive a PfVarX statistical decision.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pysam

from pfvarx.calling.evidence import SiteEvidence
from pfvarx.calling.policy import (
    DECISION_ALT_EXPLORATORY,
    DECISION_NO_CALL,
    DECISION_REFERENCE_SUPPORTED,
    SiteCall,
    minority_bin,
)
from pfvarx.config.schema import PfVarXConfig


@dataclass
class ImportResult:
    sample_id: str
    calls: list[SiteCall] = field(default_factory=list)
    n_records: int = 0
    n_with_counts: int = 0
    notes: list[str] = field(default_factory=list)


def _supported_alts(rec) -> list[tuple[int, str]]:
    """Return supported ALT alleles with their original VCF allele indexes."""
    return [
        (index, alt)
        for index, alt in enumerate(rec.alts or ())
        if alt not in (None, ".", "<NON_REF>", "*")
    ]


def _safe_get(container, key):
    """Return ``container[key]`` or ``None`` without raising on absent headers.

    pysam raises ``ValueError('Invalid header')`` when an undeclared INFO/FORMAT
    field is accessed, so presence is checked before retrieval.
    """
    try:
        if key not in container:
            return None
        return container[key]
    except (KeyError, ValueError):
        return None


def _counts(
    rec, sample: str, alt_index: int, n_alts: int
) -> tuple[int | None, int | None, int | None]:
    """Return counts for one ALT (using its zero-based ALT index), if recoverable."""
    s = rec.samples.get(sample)
    if s is None:
        return None, None, None
    ad = _safe_get(s, "AD")
    if ad is not None and len(ad) > alt_index + 1:
        ref = int(ad[0])
        alt = int(ad[alt_index + 1])
        return ref, alt, ref + alt
    ao = _safe_get(s, "AO")
    ro = _safe_get(s, "RO")
    if ao is not None and ro is not None:
        if isinstance(ao, (tuple, list)):
            if len(ao) <= alt_index:
                ao = None
            else:
                alt = int(ao[alt_index])
        elif n_alts == 1:
            alt = int(ao)
        else:
            # A scalar AO on a multiallelic record cannot be assigned to one ALT.
            ao = None
        if ao is not None:
            ref = int(ro)
            return ref, alt, ref + alt
    dp = _safe_get(s, "DP")
    if dp is None:
        dp = _safe_get(rec.info, "DP")
    if dp is not None:
        return None, None, int(dp)
    # Round-trip counts are only allele-specific for a single-ALT record.
    pvx_ref = _safe_get(rec.info, "PVX_REFCOUNT")
    pvx_alt = _safe_get(rec.info, "PVX_ALTCOUNT")
    if n_alts == 1 and pvx_ref is not None and pvx_alt is not None:
        return int(pvx_ref), int(pvx_alt), int(pvx_ref) + int(pvx_alt)
    return None, None, None


def import_variants(
    vcf_path: str,
    cfg: PfVarXConfig,
    *,
    sample_id: str,
    max_records: int | None = None,
) -> ImportResult:
    """Convert a VCF/BCF/gVCF into PfVarX SiteCall objects (no realignment)."""
    result = ImportResult(sample_id=sample_id)
    path = str(vcf_path)
    if path.endswith(".bcf"):
        opener = pysam.VariantFile(path, "rb")
    elif path.endswith(".gz"):
        opener = pysam.VariantFile(path, "r")
    else:
        opener = pysam.VariantFile(path, "r")

    with opener as vf:
        samples = list(vf.header.samples)
        if sample_id not in samples:
            if len(samples) == 1:
                result.notes.append(
                    f"Sample '{sample_id}' not in VCF; using sole sample '{samples[0]}'"
                )
                vcf_sample = samples[0]
            elif len(samples) == 0:
                result.notes.append(
                    "VCF has no sample columns; per-sample counts unavailable"
                )
                vcf_sample = sample_id
            else:
                result.notes.append(
                    f"Sample '{sample_id}' not in VCF and multiple samples present; "
                    "counts cannot be attributed"
                )
                vcf_sample = sample_id
        else:
            vcf_sample = sample_id

        for rec in vf:
            if max_records is not None and result.n_records >= max_records:
                break
            result.n_records += 1
            alts = _supported_alts(rec)
            if not alts:
                continue
            sample_record = rec.samples.get(vcf_sample)
            ao_values = _safe_get(sample_record, "AO") if sample_record is not None else None
            ro_value = _safe_get(sample_record, "RO") if sample_record is not None else None
            ad_values = _safe_get(sample_record, "AD") if sample_record is not None else None
            if (
                len(rec.alts or ()) > 1
                and ad_values is None
                and ao_values is not None
                and not isinstance(ao_values, (tuple, list))
                and ro_value is not None
            ):
                note = (
                    "MULTIALLELIC_AO_AMBIGUOUS: scalar AO cannot be assigned to individual "
                    "ALT alleles; imported without ALT counts"
                )
                if note not in result.notes:
                    result.notes.append(note)
            ref = rec.ref
            for alt_index, alt in alts:
                ref_count, alt_count, depth = _counts(
                    rec, vcf_sample, alt_index, len(rec.alts or ())
                )
                if ref_count is not None and alt_count is not None:
                    result.n_with_counts += 1
                    informative = ref_count + alt_count
                    af = alt_count / informative if informative else 0.0
                else:
                    informative = depth or 0
                    af = None

                ev = SiteEvidence(
                    sample_id=sample_id,
                    chrom=rec.chrom,
                    pos=rec.pos,
                    ref=ref,
                    alt=alt,
                    var_type="SNP" if len(ref) == len(alt) == 1 else "INDEL",
                    ref_count=ref_count or 0,
                    alt_count=alt_count or 0,
                    other_count=0,
                    depth=depth or 0,
                    informative_depth=informative,
                    alt_fraction=af,
                    sequencing_mode=cfg.run.mode,
                    umi_status="NOT_APPLICABLE" if cfg.run.mode != "tbs_umi" else "IMPORTED",
                    notes=["IMPORTED_CALLSET"],
                )
                if af is None:
                    decision = DECISION_NO_CALL
                    filters = ["UncalibratedModel"]
                elif alt_count == 0:
                    decision = DECISION_REFERENCE_SUPPORTED
                    filters = []
                else:
                    # Imported calls are exploratory by definition: PfVarX did not
                    # compute their likelihoods, so they can never be PASS.
                    decision = DECISION_ALT_EXPLORATORY
                    filters = ["UncalibratedModel"]
                call = SiteCall(
                    evidence=ev,
                    decision=decision,
                    filters=filters,
                    bin=minority_bin(af, cfg) if af is not None else "BELOW_THRESHOLD",
                    model="imported",
                    model_version="import-v1",
                    calibration_profile="IMPORTED",
                    policy_status=cfg.run.policy_status,
                    policy_id=cfg.policy_id(),
                    notes=["Imported from external callset; PfVarX likelihoods not computed"],
                )
                result.calls.append(call)

    if result.n_with_counts == 0 and result.n_records:
        result.notes.append("IMPORTED_COUNTS_UNAVAILABLE: no unambiguous AD/AO/RO counts found")
    return result
