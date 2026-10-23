"""Machine-readable release status registry.

The specification (section 17) requires explicit, non-inflated status labels.
A passing test suite is *not* evidence of biological accuracy; the two axes are
tracked separately here.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Literal

Status = Literal[
    "NOT_IMPLEMENTED",
    "IMPLEMENTED_NOT_VALIDATED",
    "TESTED_NOT_SCIENTIFICALLY_VALIDATED",
    "VALIDATED_WITH_RESTRICTIONS",
    "RESEARCH_USE_READY",
]


@dataclass(frozen=True)
class ComponentStatus:
    component: str
    software_status: Status
    scientific_status: Status
    note: str


_COMPONENTS: tuple[ComponentStatus, ...] = (
    ComponentStatus(
        "configuration",
        "TESTED_NOT_SCIENTIFICALLY_VALIDATED",
        "IMPLEMENTED_NOT_VALIDATED",
        "Pydantic schema, bounds, profile precedence, config hashing, input-only initialization, "
        "and configurable results/log/work roots.",
    ),
    ComponentStatus(
        "input_validation_routing",
        "TESTED_NOT_SCIENTIFICALLY_VALIDATED",
        "IMPLEMENTED_NOT_VALIDATED",
        "FASTQ/BAM/CRAM/VCF/BCF/gVCF/FASTA detection and route selection.",
    ),
    ComponentStatus(
        "preprocessing_alignment",
        "TESTED_NOT_SCIENTIFICALLY_VALIDATED",
        "NOT_IMPLEMENTED",
        "Thin, versioned wrappers around fastp/BWA-MEM2/samtools/Picard/fgbio; "
        "synthetic paired-FASTQ workflow smoke-tested; external binaries are not "
        "bundled in the Python wheel and real-data accuracy is not validated.",
    ),
    ComponentStatus(
        "evidence_extraction",
        "TESTED_NOT_SCIENTIFICALLY_VALIDATED",
        "IMPLEMENTED_NOT_VALIDATED",
        "CIGAR-aware pileup, fragment de-duplication, strand/position/context metrics.",
    ),
    ComponentStatus(
        "statistical_models",
        "TESTED_NOT_SCIENTIFICALLY_VALIDATED",
        "IMPLEMENTED_NOT_VALIDATED",
        "Composite-AF binomial/beta-binomial models and exact one-sided count-tail "
        "p-values with complete-family Bonferroni adjustment; null calibration "
        "requires user-supplied controls.",
    ),
    ComponentStatus(
        "decision_policy",
        "TESTED_NOT_SCIENTIFICALLY_VALIDATED",
        "IMPLEMENTED_NOT_VALIDATED",
        "Provisional thresholds; validated PASS fails closed without verifiable "
        "assay-validation/LoD provenance.",
    ),
    ComponentStatus(
        "annotation",
        "TESTED_NOT_SCIENTIFICALLY_VALIDATED",
        "NOT_IMPLEMENTED",
        "Internal GFF3/CDS consequence engine and conservative marker matching; "
        "SnpEff is not integrated into per-call annotation.",
    ),
    ComponentStatus(
        "downstream_analyses",
        "IMPLEMENTED_NOT_VALIDATED",
        "NOT_IMPLEMENTED",
        "Diversity, mixed-infection summaries implemented; DEploid/hmmIBD/IQ-TREE are "
        "guarded wrappers that fail closed when assumptions or binaries are absent.",
    ),
    ComponentStatus(
        "unified_qc",
        "TESTED_NOT_SCIENTIFICALLY_VALIDATED",
        "IMPLEMENTED_NOT_VALIDATED",
        "Metric ontology with value/unit/source/status/reason/policy_id.",
    ),
    ComponentStatus(
        "dashboard_reporting",
        "IMPLEMENTED_NOT_VALIDATED",
        "NOT_IMPLEMENTED",
        "File-backed Streamlit dashboard and CSV/XLSX/PNG/PDF exporters.",
    ),
    ComponentStatus(
        "benchmarking",
        "TESTED_NOT_SCIENTIFICALLY_VALIDATED",
        "NOT_IMPLEMENTED",
        "Simulation + ablation framework executes; no external truth set is bundled.",
    ),
)


RELEASE_STATUS: dict[str, dict[str, str]] = {c.component: asdict(c) for c in _COMPONENTS}


def get_status(component: str | None = None) -> dict:
    """Return the release status registry, optionally for one component."""
    if component is None:
        return RELEASE_STATUS
    if component not in RELEASE_STATUS:
        raise KeyError(f"Unknown component {component!r}")
    return RELEASE_STATUS[component]


def overall_status() -> Status:
    """Conservative aggregate: the weakest scientific status across components."""
    order: list[Status] = [
        "NOT_IMPLEMENTED",
        "IMPLEMENTED_NOT_VALIDATED",
        "TESTED_NOT_SCIENTIFICALLY_VALIDATED",
        "VALIDATED_WITH_RESTRICTIONS",
        "RESEARCH_USE_READY",
    ]
    worst = min(order.index(c.scientific_status) for c in _COMPONENTS)
    return order[worst]


def status_banner() -> str:
    return (
        "RESEARCH USE ONLY — validation status: "
        f"{overall_status()} (no component is clinically validated)"
    )
