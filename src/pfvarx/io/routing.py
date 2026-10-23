"""Configuration-driven input routing.

Routing is explicit and auditable: a valid existing BAM/CRAM is never realigned
merely because the workflow started, and a VCF/BCF/gVCF follows the variant
import route rather than read-level calling.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

from pfvarx.config.schema import PfVarXConfig
from pfvarx.errors import RoutingError
from pfvarx.io.validators import SampleRecord, detect_type

ROUTE_FASTQ = "fastq_alignment"
ROUTE_BAM = "alignment_reuse"
ROUTE_VARIANT = "variant_import"

STAGE_ORDER = [
    "validate",
    "preprocess",
    "align",
    "bam_qc",
    "discover",
    "evidence",
    "call",
    "normalize",
    "variant_qc",
    "annotate",
    "downstream",
    "unified_qc",
    "report",
]


@dataclass
class RouteDecision:
    sample_id: str
    route: str
    input_type: str
    evidence_path: str
    stages: list[str] = field(default_factory=list)
    reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)


def _stages_for_route(route: str, cfg: PfVarXConfig) -> list[str]:
    if route == ROUTE_FASTQ:
        stages = ["validate", "preprocess", "align", "bam_qc", "discover", "evidence", "call",
                  "normalize", "variant_qc"]
    elif route == ROUTE_BAM:
        stages = ["validate", "bam_qc", "discover", "evidence", "call", "normalize", "variant_qc"]
    elif route == ROUTE_VARIANT:
        stages = ["validate", "normalize", "variant_qc"]
    else:  # pragma: no cover - defensive
        raise RoutingError(f"Unknown route {route}")
    if cfg.annotation.enabled:
        stages.append("annotate")
    if cfg.downstream.enabled:
        stages.append("downstream")
    stages.extend(["unified_qc", "report"])
    return stages


def route_sample(sample: SampleRecord, cfg: PfVarXConfig) -> RouteDecision:
    declared = sample.input_type
    evidence = sample.evidence_path()
    if evidence is None:
        raise RoutingError(
            f"Sample {sample.sample_id} declares input_type={declared} but no evidence path",
            sample=sample.sample_id,
            action="Populate fastq_1, bam or vcf in the samplesheet",
        )
    if not Path(evidence).exists():
        raise RoutingError(
            f"Sample {sample.sample_id}: input file not found: {evidence}",
            sample=sample.sample_id,
            action="Correct the path in the samplesheet (fastq_1/fastq_2, bam or vcf).",
        )
    detected = detect_type(evidence)
    warnings: list[str] = []
    reasons: list[str] = []

    if (
        declared != "auto"
        and declared != detected
        and {declared, detected} not in ({"bam", "cram"},)
    ):
        raise RoutingError(
            f"Sample {sample.sample_id}: declared input_type={declared} but detected {detected}",
            sample=sample.sample_id,
            cause="File extension/content disagree",
            action="Correct input_type or supply the intended file",
        )

    if detected == "fastq":
        if not sample.fastq_1:
            raise RoutingError(f"Sample {sample.sample_id}: FASTQ route requires fastq_1", sample=sample.sample_id)
        route = ROUTE_FASTQ
        reasons.append("FASTQ evidence requires trimming and alignment with BWA-MEM2")
        if sample.fastq_2:
            reasons.append("Paired-end layout detected")
        else:
            warnings.append("Single-end FASTQ: fragment-level de-duplication is limited")
    elif detected in {"bam", "cram"}:
        route = ROUTE_BAM
        reasons.append(f"Existing {detected.upper()} reused; no realignment")
        if not cfg.inputs.accept_existing_bam:
            raise RoutingError(
                f"Sample {sample.sample_id}: existing BAM/CRAM supplied but accept_existing_bam=false",
                sample=sample.sample_id,
                action="Set inputs.accept_existing_bam=true or supply FASTQ",
            )
    elif detected in {"vcf", "bcf", "gvcf"}:
        route = ROUTE_VARIANT
        reasons.append(f"{detected.upper()} imported through the variant-input route")
        warnings.append(
            "Imported variants are not equivalent to newly called PfVarX evidence; "
            "provenance and quality are assessed, not re-derived."
        )
        if detected == "gvcf":
            warnings.append("gVCF: symbolic NON_REF blocks are not treated as alternate alleles")
    else:
        raise RoutingError(
            f"Sample {sample.sample_id}: unsupported input type {detected!r}",
            sample=sample.sample_id,
            action="Supply FASTQ, BAM, CRAM, VCF, BCF or gVCF",
        )

    if cfg.run.mode == "tbs_umi" and route == ROUTE_FASTQ:
        warnings.append("UMI mode: UMI extraction/consensus must be validated for this assay before molecule-level evidence is used")
    if cfg.run.mode == "tbs_no_umi" and route == ROUTE_BAM:
        warnings.append("TBS without UMI: duplicate marking is intentionally not applied")

    return RouteDecision(
        sample_id=sample.sample_id,
        route=route,
        input_type=detected,
        evidence_path=str(evidence),
        stages=_stages_for_route(route, cfg),
        reasons=reasons,
        warnings=warnings,
    )


def plan_run(samples: list[SampleRecord], cfg: PfVarXConfig) -> list[RouteDecision]:
    if not samples:
        raise RoutingError("No samples to route", action="Populate the samplesheet")
    decisions = [route_sample(s, cfg) for s in samples]
    ids = [d.sample_id for d in decisions]
    if len(set(ids)) != len(ids):
        raise RoutingError("Duplicate sample IDs after routing")
    return decisions


def route_all(
    samples: list[SampleRecord], cfg: PfVarXConfig
) -> tuple[list[RouteDecision], list[dict]]:
    """Route every sample, collecting failures instead of stopping at the first.

    Returns ``(decisions, errors)``.  ``errors`` holds one structured record per
    sample that could not be routed, so a validation report can name *every*
    problem in one pass rather than making the user fix them one at a time.
    """
    if not samples:
        raise RoutingError("No samples to route", action="Populate the samplesheet")

    decisions: list[RouteDecision] = []
    errors: list[dict] = []
    for sample in samples:
        try:
            decisions.append(route_sample(sample, cfg))
        except RoutingError as exc:
            payload = exc.as_dict()
            payload.setdefault("sample", sample.sample_id)
            errors.append(payload)

    ids = [d.sample_id for d in decisions]
    if len(set(ids)) != len(ids):
        errors.append(
            {
                "error_type": "RoutingError",
                "stage": "routing",
                "message": "Duplicate sample IDs after routing",
                "corrective_action": "Give every row a unique sample_id.",
            }
        )
    return decisions, errors


def summarise_plan(decisions: list[RouteDecision]) -> dict:
    from collections import Counter

    routes = Counter(d.route for d in decisions)
    return {
        "n_samples": len(decisions),
        "routes": dict(routes),
        "samples": [d.as_dict() for d in decisions],
        "requires_alignment": routes.get(ROUTE_FASTQ, 0),
        "reuses_alignment": routes.get(ROUTE_BAM, 0),
        "imports_variants": routes.get(ROUTE_VARIANT, 0),
    }
