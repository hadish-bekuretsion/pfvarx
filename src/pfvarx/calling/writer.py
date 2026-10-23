"""VCF/BCF and structured-evidence writers.

Output contract (spec section 7):

* The VCF/BCF carries the variant record and a compact INFO summary.
* A **separate** structured evidence table carries the full per-site evidence:
  reference/alternate/uninformative counts, depth, allele fraction, model and
  version, log-likelihoods, LRT statistic and p/q values, calibrated posterior
  (if any), calibration profile, filter status and reasons, uncertainty flags,
  sequencing mode and UMI status.
* ``QUAL`` is emitted as ``.`` unless a genuinely calibrated posterior exists.
  A raw likelihood-ratio statistic is never placed in ``QUAL``.
* Every record carries the policy id so that a call can be traced to the exact
  configuration that produced it.
"""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path

import pysam

from pfvarx.calling.policy import DECISION_ALT_PASS, POLICY_VERSION, SiteCall
from pfvarx.config.schema import PfVarXConfig
from pfvarx.errors import OutputError

EVIDENCE_COLUMNS = [
    "sample_id", "chrom", "pos", "ref", "alt", "var_type",
    "ref_count", "alt_count", "other_count", "depth", "informative_depth",
    "alt_fraction", "af_ci_low", "af_ci_high",
    "mean_base_quality_alt", "mean_base_quality_ref", "mean_mapping_quality",
    "alt_forward", "alt_reverse", "ref_forward", "ref_reverse",
    "strand_bias_pvalue", "homopolymer_length", "left_aligned",
    "duplicate_fragments_excluded", "low_mapq_excluded", "low_bq_excluded",
    "secondary_excluded", "sequencing_mode", "umi_status",
    "decision", "filter_status", "filter_reasons", "minority_bin",
    "model", "model_version", "log_likelihood_h0", "log_likelihood_h1",
    "lrt_statistic", "lrt_pvalue", "lrt_pvalue_adjusted", "posterior_h1",
    "calibration_profile", "policy_status", "policy_id", "policy_notes",
]

_FILTERS = (
    "LowEffectiveDepth", "LowAltSupport", "UncalibratedModel",
    "StrandBiasReview", "HomopolymerIndelReview", "OutsideCallableRegion",
    "AmbiguousMapping", "UMIModelInvalid", "ReferenceMismatch",
)


def _vcf_header(cfg: PfVarXConfig, sample_id: str) -> pysam.VariantHeader:
    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.3")
    header.add_meta("source", value=f"PfVarX {POLICY_VERSION}")
    header.add_meta("pfvarx_policy_id", value=cfg.policy_id())
    header.add_meta("pfvarx_policy_status", value=cfg.run.policy_status)
    header.add_meta("pfvarx_sequencing_mode", value=cfg.run.mode)
    header.add_meta(
        "pfvarx_status_banner",
        value="RESEARCH USE ONLY - see release status in the evidence table",
    )
    info = [
        ("PVX_DECISION", "1", "String", "PfVarX decision category"),
        ("PVX_BIN", "1", "String", "Minority fraction bin"),
        ("PVX_MODEL", "1", "String", "Count model used"),
        ("PVX_MODEL_VERSION", "1", "String", "Model/policy version"),
        ("PVX_REFCOUNT", "1", "Integer", "Reference-supporting observations"),
        ("PVX_ALTCOUNT", "1", "Integer", "Alternate-supporting observations"),
        ("PVX_DEPTH", "1", "Integer", "Informative depth"),
        ("PVX_AF", "1", "Float", "Alternate allele fraction"),
        ("PVX_AF_CI", "2", "Float", "Allele-fraction credible interval"),
        ("PVX_LRT", "1", "Float", "Likelihood-ratio statistic (NOT a posterior)"),
        ("PVX_LRT_P", "1", "Float", "LRT p-value"),
        ("PVX_Q", "1", "Float", "BH-adjusted q-value"),
        ("PVX_POSTERIOR", "1", "Float", "Calibrated posterior P(H1) if available"),
        ("PVX_CALIBRATION", "1", "String", "Calibration profile"),
        ("PVX_POLICY_STATUS", "1", "String", "Policy validation status"),
        ("PVX_POLICY_ID", "1", "String", "Configuration policy identifier"),
        ("PVX_UMI_STATUS", "1", "String", "UMI consensus status"),
    ]
    for i, n, t, d in info:
        header.add_meta(
            "INFO",
            items=[("ID", i), ("Number", n), ("Type", t), ("Description", d)],
        )
    for f in _FILTERS:
        header.add_meta(
            "FILTER",
            items=[("ID", f), ("Description", f"PfVarX review filter {f}")],
        )
    header.add_meta(
        "FORMAT",
        items=[("ID", "GT"), ("Number", "1"), ("Type", "String"), ("Description", "Genotype")],
    )
    header.contigs.add(cfg.reference.assembly_id)
    header.add_sample(sample_id)
    return header


def _qual(call: SiteCall, cfg: PfVarXConfig) -> float | None:
    """Return a QUAL only when a calibrated posterior is genuinely available."""
    if (
        cfg.run.policy_status == "VALIDATED"
        and call.posterior_h1 is not None
        and call.decision == DECISION_ALT_PASS
    ):
        return round(-10.0 * math.log10(max(1e-300, 1.0 - call.posterior_h1)), 2)
    return None


def write_vcf(
    calls: list[SiteCall],
    cfg: PfVarXConfig,
    out_path: str,
    *,
    sample_id: str,
    bcf: bool = False,
) -> str:
    """Write calls to VCF or BCF.  Returns the path written."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    header = _vcf_header(cfg, sample_id)
    for chrom in sorted({c.evidence.chrom for c in calls}):
        if chrom not in header.contigs:
            header.contigs.add(chrom)
    mode = "wb" if bcf else "w"
    tmp = out.with_suffix(out.suffix + ".tmp")
    try:
        with pysam.VariantFile(str(tmp), mode, header=header) as vf:
            for call in sorted(calls, key=lambda c: (c.evidence.chrom, c.evidence.pos)):
                ev = call.evidence
                rec = vf.new_record(
                    contig=ev.chrom,
                    start=ev.pos - 1,
                    stop=ev.pos - 1 + len(ev.ref),
                    alleles=(ev.ref, ev.alt),
                    id=None,
                    qual=_qual(call, cfg),
                    filter=list(call.filters) if call.filters else "PASS",
                )
                rec.info["PVX_DECISION"] = call.decision
                rec.info["PVX_BIN"] = call.bin
                rec.info["PVX_MODEL"] = call.model
                rec.info["PVX_MODEL_VERSION"] = call.model_version
                rec.info["PVX_REFCOUNT"] = ev.ref_count
                rec.info["PVX_ALTCOUNT"] = ev.alt_count
                rec.info["PVX_DEPTH"] = ev.informative_depth
                rec.info["PVX_AF"] = round(ev.alt_fraction, 6) if ev.alt_fraction is not None else None
                if ev.af_ci_low is not None and ev.af_ci_high is not None:
                    rec.info["PVX_AF_CI"] = (round(ev.af_ci_low, 6), round(ev.af_ci_high, 6))
                if call.lrt_statistic is not None:
                    rec.info["PVX_LRT"] = round(call.lrt_statistic, 4)
                if call.lrt_pvalue is not None:
                    rec.info["PVX_LRT_P"] = round(call.lrt_pvalue, 8)
                if call.lrt_pvalue_adjusted is not None:
                    rec.info["PVX_Q"] = round(call.lrt_pvalue_adjusted, 8)
                if call.posterior_h1 is not None:
                    rec.info["PVX_POSTERIOR"] = round(call.posterior_h1, 6)
                rec.info["PVX_CALIBRATION"] = call.calibration_profile
                rec.info["PVX_POLICY_STATUS"] = call.policy_status
                rec.info["PVX_POLICY_ID"] = call.policy_id or cfg.policy_id()
                rec.info["PVX_UMI_STATUS"] = ev.umi_status
                # P. falciparum is haploid in the intra-erythrocytic stage, so
                # genotypes are single-allele (haploid) calls, not diploid.
                rec.samples[sample_id]["GT"] = (1,) if call.is_alt else (0,)
                vf.write(rec)
        os.replace(tmp, out)
    except Exception as exc:  # pragma: no cover - surfaced as domain error
        if tmp.exists():
            tmp.unlink(missing_ok=True)
        raise OutputError(
            "Failed to write variant output",
            cause=str(exc),
            action="Check output directory permissions and reference contig names.",
        ) from exc
    return str(out)


def write_evidence_table(calls: list[SiteCall], out_path: str) -> str:
    """Write the full structured evidence table as TSV."""
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=EVIDENCE_COLUMNS, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for call in sorted(calls, key=lambda c: (c.evidence.chrom, c.evidence.pos)):
            row = call.as_dict()
            writer.writerow({k: row.get(k, "") for k in EVIDENCE_COLUMNS})
    os.replace(tmp, out)
    return str(out)
