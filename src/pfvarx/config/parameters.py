"""Quantitative parameter registry.

Every tunable in PfVarX is declared here with its unit, datatype, accepted
range, provisional default, sequencing-mode applicability, rationale, source,
policy status, calibration-dataset identifier and effective version.

Policy status vocabulary (specification section 20):

``PROVISIONAL_DEVELOPMENT``
    Engineering starting point. Not a validated cutoff. Results derived from it
    must be labelled ``UNVALIDATED``.
``VALIDATED``
    Backed by a signed, reproducible, mode-specific validation report. Requires
    ``calibration_dataset_id``.
``DERIVED``
    Computed at runtime from data (e.g. median depth); never user-set.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

MODES = ("wgs", "tbs_no_umi", "tbs_umi", "all")

POLICY_PROVISIONAL = "PROVISIONAL_DEVELOPMENT"
POLICY_VALIDATED = "VALIDATED"
POLICY_DERIVED = "DERIVED"


@dataclass(frozen=True)
class Parameter:
    name: str
    unit: str
    datatype: str
    default: Any
    minimum: float | None = None
    maximum: float | None = None
    allowed: tuple[Any, ...] | None = None
    modes: tuple[str, ...] = ("all",)
    rationale: str = ""
    source: str = ""
    policy_status: str = POLICY_PROVISIONAL
    calibration_dataset_id: str | None = None
    version: str = "1.0.0"
    category: str = "general"

    def coerce(self, value: Any) -> Any:
        if self.datatype == "float":
            return float(value)
        if self.datatype == "int":
            if isinstance(value, bool):
                raise TypeError(f"{self.name}: boolean is not a valid integer")
            return int(value)
        if self.datatype == "bool":
            if isinstance(value, bool):
                return value
            if isinstance(value, str) and value.lower() in {"true", "false"}:
                return value.lower() == "true"
            raise TypeError(f"{self.name}: expected boolean, got {value!r}")
        if self.datatype == "str":
            return str(value)
        return value

    def check(self, value: Any) -> Any:
        value = self.coerce(value)
        if self.allowed is not None and value not in self.allowed:
            raise ValueError(f"{self.name}: {value!r} not in allowed values {self.allowed}")
        if self.datatype in {"float", "int"}:
            if self.minimum is not None and value < self.minimum:
                raise ValueError(f"{self.name}: {value} < minimum {self.minimum}")
            if self.maximum is not None and value > self.maximum:
                raise ValueError(f"{self.name}: {value} > maximum {self.maximum}")
        return value

    def applies_to(self, mode: str) -> bool:
        return "all" in self.modes or mode in self.modes

    def as_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# Registry. Sources are real, verifiable publications / tool documentation.
# ---------------------------------------------------------------------------
_REGISTRY: list[Parameter] = [
    # -- reference and inputs ------------------------------------------------
    Parameter(
        "reference.assembly_id", "identifier", "str", "REQUIRED_PINNED_ASSEMBLY_ID",
        modes=("all",), category="reference",
        rationale="Assembly identity must be pinned; contig names alone do not identify a build.",
        source="MalariaGEN Pf7 bioinformatics methods (2023)",
    ),
    Parameter(
        "reference.sha256", "hex digest", "str", "REQUIRED_SHA256", category="reference",
        rationale="Detects silent reference substitution and CRAM/BAM mismatch.",
        source="SAMtools/CRAM specification; GA4GH refget",
    ),
    Parameter(
        "reference.nuclear_ploidy", "copies", "int", 1, minimum=1, maximum=1, category="reference",
        rationale="P. falciparum is haploid in the human blood stage; allele fractions must not be forced to 0/1 in polyclonal samples.",
        source="Gardner et al. (2002) Nature 419:498-511",
    ),
    Parameter(
        "preprocessing.fastp_min_phred", "Phred", "int", 20, minimum=0, maximum=40,
        modes=("all",), category="preprocessing",
        rationale="Q20 = nominal 1% base-call error; a trimming heuristic, not an empirical error rate.",
        source="fastp documentation (Chen et al. 2018, Bioinformatics 34:884-890)",
    ),
    Parameter(
        "preprocessing.fastp_min_length", "bp", "int", 50, minimum=0, maximum=1000,
        modes=("all",), category="preprocessing",
        rationale="Reads shorter than 50 bp map ambiguously in AT-rich P. falciparum.",
        source="MalariaGEN Pf7 methods; BWA-MEM2 recommendations",
    ),
    Parameter(
        "preprocessing.fastp_max_n_fraction", "fraction", "float", 0.05, minimum=0.0, maximum=1.0,
        modes=("all",), category="preprocessing",
        rationale="Excess N calls indicate quality collapse.",
        source="fastp documentation",
    ),
    Parameter(
        "preprocessing.mapping_quality_min", "MAPQ", "int", 30, minimum=0, maximum=60,
        modes=("all",), category="preprocessing",
        rationale="MAPQ is not a calibrated probability of correct mapping; 30 is a conservative read-level filter.",
        source="Li (2013) arXiv:1303.3997; MalariaGEN Pf7 methods",
    ),
    Parameter(
        "preprocessing.base_quality_min", "Phred", "int", 20, minimum=0, maximum=40,
        modes=("all",), category="preprocessing",
        rationale="Excludes low-confidence bases from allele counting; sensitivity at Q30 should be reported.",
        source="MalariaGEN Pf7 methods",
    ),
    Parameter(
        "preprocessing.count_overlapping_mates_once", "boolean", "bool", True,
        modes=("all",), category="preprocessing",
        rationale="Overlapping mates are one physical fragment, not two independent observations.",
        source="samtools mpileup -x / GATK read-overlap conventions",
    ),
    Parameter(
        "preprocessing.wgs_mark_duplicates", "boolean", "bool", True,
        modes=("wgs",), category="preprocessing",
        rationale="PCR duplicates inflate depth in WGS and violate independence assumptions.",
        source="Picard MarkDuplicates documentation",
    ),
    Parameter(
        "preprocessing.duplicate_exclusion_from_evidence", "boolean", "bool", True,
        modes=("wgs",), category="preprocessing",
        rationale="Marked duplicates are excluded from independent evidence but retained on disk.",
        source="MalariaGEN Pf7 methods",
    ),
    # -- callability ---------------------------------------------------------
    Parameter(
        "callability.wgs_min_independent_depth", "fragments", "int", 10, minimum=1, maximum=100000,
        modes=("wgs",), category="callability",
        rationale="Provisional majority-call evaluation depth; low depth is not proof of absence.",
        source="MalariaGEN Pf7 methods; Pf7 QC description (2023)",
    ),
    Parameter(
        "callability.targeted_min_independent_depth", "fragments", "int", 100, minimum=1, maximum=1000000,
        modes=("tbs_no_umi",), category="callability",
        rationale="Exploratory minority-allele screening depth; 100x total depth does not establish 1% LoD.",
        source="Reporting considerations for targeted malaria molecular surveillance (2025)",
    ),
    Parameter(
        "callability.umi_min_consensus_molecules_majority", "molecules", "int", 20,
        minimum=1, maximum=100000, modes=("tbs_umi",), category="callability",
        rationale="Molecule-level support for majority-allele assessment.",
        source="fgbio documentation; UMI consensus literature",
    ),
    Parameter(
        "callability.umi_min_consensus_molecules_exploratory", "molecules", "int", 100,
        minimum=1, maximum=1000000, modes=("tbs_umi",), category="callability",
        rationale="Exploratory minority screening at molecule level.",
        source="fgbio documentation",
    ),
    Parameter(
        "callability.wgs_high_depth_multiplier_for_review", "multiplier", "float", 5.0,
        minimum=1.0, maximum=100.0, modes=("wgs",), category="callability",
        rationale="Depth above max(500, 5x median) flags possible mapping collapse; a review flag, not a fail.",
        source="MalariaGEN Pf7 methods",
    ),
    Parameter(
        "callability.wgs_max_depth_floor", "fragments", "int", 500, minimum=0, maximum=100000,
        modes=("wgs",), category="callability",
        rationale="Absolute floor of the high-depth review rule.",
        source="MalariaGEN Pf7 methods",
    ),
    Parameter(
        "callability.wgs_min_callable_fraction", "fraction", "float", 0.80, minimum=0.0, maximum=1.0,
        modes=("wgs",), category="callability",
        rationale="Sample-level warning when <80% of the callable mask reaches 10x.",
        source="MalariaGEN Pf7 methods",
    ),
    # -- candidate discovery -------------------------------------------------
    Parameter(
        "candidate.snp_min_alt_observations", "fragments", "int", 2, minimum=1, maximum=100,
        modes=("all",), category="candidate",
        rationale="Permissive discovery only; two supporting fragments are not report-level evidence.",
        source="Development policy; cf. BCFtools mpileup default",
    ),
    Parameter(
        "candidate.indel_min_alt_observations", "fragments", "int", 3, minimum=1, maximum=100,
        modes=("all",), category="candidate",
        rationale="INDELs require more support because alignment ambiguity inflates false positives.",
        source="Development policy; cf. FreeBayes/GATK indel defaults",
    ),
    Parameter(
        "candidate.exploratory_min_alt_fraction", "fraction", "float", 0.01, minimum=0.0, maximum=1.0,
        modes=("all",), category="candidate",
        rationale="Discovery floor; below-threshold evidence is recorded, not silently dropped.",
        source="Development policy",
    ),
    Parameter(
        "candidate.max_small_indel_bp", "bp", "int", 50, minimum=1, maximum=50,
        modes=("all",), category="candidate",
        rationale="Scope limit: SNPs and small INDELs only.",
        source="Specification section 2",
    ),
    Parameter(
        "candidate.homopolymer_review_length", "bp", "int", 6, minimum=1, maximum=50,
        modes=("all",), category="candidate",
        rationale="Homopolymer run >=6 bp is a review flag requiring INDEL-specific calibration.",
        source="Development policy; homopolymer indel literature",
    ),
    Parameter(
        "candidate.strand_bias_alpha", "probability", "float", 0.01, minimum=0.0, maximum=1.0,
        modes=("all",), category="candidate",
        rationale="Binomial strand test; unreliable with sparse data, so used only as a review flag.",
        source="Guo et al. (2012) strand-bias assessment",
    ),
    # -- statistics ----------------------------------------------------------
    Parameter(
        "statistics.count_model", "enum", "str", "auto", allowed=("binomial", "beta_binomial", "auto"),
        modes=("all",), category="statistics",
        rationale="Model selection is data-driven; beta-binomial only when calibration demonstrates overdispersion.",
        source="Development policy; beta-binomial overdispersion literature",
    ),
    Parameter(
        "statistics.error_rate", "probability", "float", None, minimum=0.0, maximum=1.0,
        modes=("all",), category="statistics",
        rationale="Empirical per-mode error probability. Without calibration, sensitivity scenarios 0.001/0.005/0.01 are used and results marked unvalidated.",
        source="Development policy; requires assay-specific control data",
    ),
    Parameter(
        "statistics.beta_binomial_rho", "fraction", "float", None, minimum=0.0, maximum=0.999999,
        modes=("all",), category="statistics",
        rationale="Overdispersion; rho=0 is the binomial limit. Must be fitted from independent calibration data.",
        source="Development policy",
    ),
    Parameter(
        "statistics.min_alt_fraction_h1", "fraction", "float", 0.01, minimum=0.0, maximum=1.0,
        modes=("all",), category="statistics",
        rationale="Biological detection boundary f_min used to define H1.",
        source="Development policy; must be replaced by validated LoD",
    ),
    Parameter(
        "statistics.lrt_null_calibration", "enum", "str", "bootstrap", allowed=("bootstrap", "chi2_mixture", "none"),
        modes=("all",), category="statistics",
        rationale="The null lies on the boundary and sites are selected, so chi-square asymptotics are not assumed by default.",
        source="Self & Liang (1987); parametric bootstrap literature",
    ),
    Parameter(
        "statistics.bootstrap_replicates", "replicates", "int", 200, minimum=0, maximum=100000,
        modes=("all",), category="statistics",
        rationale="Null calibration replicates; 200 is a development default for speed.",
        source="Development policy",
    ),
    Parameter(
        "statistics.exploratory_fdr_q", "q-value", "float", 0.05, minimum=0.0, maximum=1.0,
        modes=("all",), category="statistics",
        rationale="Benjamini-Hochberg exploratory FDR target; a q-value is not a posterior probability.",
        source="Benjamini & Hochberg (1995) JRSS-B 57:289-300",
    ),
    Parameter(
        "statistics.sensitivity_fdr_q", "q-value", "float", 0.01, minimum=0.0, maximum=1.0,
        modes=("all",), category="statistics",
        rationale="Stringent sensitivity setting reported alongside the primary FDR target.",
        source="Benjamini & Hochberg (1995)",
    ),
    Parameter(
        "statistics.bayesian_enabled", "boolean", "bool", False,
        modes=("all",), category="statistics",
        rationale="Priors are disabled by default; enabling requires documented, defensible priors.",
        source="Specification section 6.2 Model D",
    ),
    Parameter(
        "statistics.bayesian_prior_h1", "probability", "float", 0.01, minimum=0.0, maximum=1.0,
        modes=("all",), category="statistics",
        rationale="Prior probability of a true alternate allele; must be justified, never used to force known alleles.",
        source="Development policy",
    ),
    Parameter(
        "statistics.posterior_threshold_if_validated", "probability", "float", None,
        minimum=0.0, maximum=1.0, modes=("all",), category="statistics",
        rationale="Only usable as a candidate gate after calibration; cannot replace negative controls.",
        source="Development policy",
    ),
    Parameter(
        "statistics.low_depth_fallback", "enum", "str", "conservative",
        allowed=("conservative", "fail"), modes=("all",), category="statistics",
        rationale="When calibration is insufficient, fail closed for validated calls but emit labelled exploratory output.",
        source="Specification section 20.5",
    ),
    # -- reporting / minority bins ------------------------------------------
    Parameter(
        "reporting.minority_candidate_fraction", "fraction", "float", 0.05, minimum=0.0, maximum=1.0,
        modes=("all",), category="reporting",
        rationale="Lower bound of MINORITY_CANDIDATE descriptive bin; not a clone count or validated LoD.",
        source="Development policy",
    ),
    Parameter(
        "reporting.intermediate_fraction", "fraction", "float", 0.20, minimum=0.0, maximum=1.0,
        modes=("all",), category="reporting",
        rationale="Lower bound of INTERMEDIATE_FRACTION bin.",
        source="Development policy",
    ),
    Parameter(
        "reporting.predominant_fraction", "fraction", "float", 0.80, minimum=0.0, maximum=1.0,
        modes=("all",), category="reporting",
        rationale="Lower bound of PREDOMINANT_ALT bin.",
        source="Development policy",
    ),
    Parameter(
        "reporting.pass_requires_assay_validated_lod", "boolean", "bool", True,
        modes=("all",), category="reporting",
        rationale="ALT_PASS_VALIDATED is impossible without a validated assay LoD.",
        source="Specification section 20.5",
    ),
    Parameter(
        "reporting.uncertainty_always_visible", "boolean", "bool", True,
        modes=("all",), category="reporting",
        rationale="Uncertainty flags must accompany every call in reports.",
        source="Specification section 20.8",
    ),
    # -- population / downstream --------------------------------------------
    Parameter(
        "population.sample_missingness_max", "fraction", "float", 0.20, minimum=0.0, maximum=1.0,
        modes=("all",), category="population",
        rationale="Exploratory cohort filter; sensitivity at 5/10/20%.",
        source="Development policy",
    ),
    Parameter(
        "population.variant_site_missingness_max", "fraction", "float", 0.10, minimum=0.0, maximum=1.0,
        modes=("all",), category="population",
        rationale="Exploratory PCA/structure filter.",
        source="Development policy",
    ),
    Parameter(
        "population.maf_min", "fraction", "float", 0.05, minimum=0.0, maximum=0.5,
        modes=("all",), category="population",
        rationale="Exploratory common-variant PCA; never applied to rare resistance-marker detection.",
        source="Development policy",
    ),
    Parameter(
        "population.ld_window_snps", "SNPs", "int", 50, minimum=1, maximum=100000,
        modes=("all",), category="population",
        rationale="Provisional LD-pruning window; only valid if haploid/mixed encoding assumptions hold.",
        source="PLINK --indep-pairwise convention",
    ),
    Parameter(
        "population.ld_step_snps", "SNPs", "int", 5, minimum=1, maximum=10000,
        modes=("all",), category="population",
        rationale="LD-pruning step size.",
        source="PLINK --indep-pairwise convention",
    ),
    Parameter(
        "population.ld_r2_max", "r^2", "float", 0.2, minimum=0.0, maximum=1.0,
        modes=("all",), category="population",
        rationale="LD-pruning r^2 threshold; disable if haploid assumptions are invalid.",
        source="PLINK --indep-pairwise convention",
    ),
    Parameter(
        "phylogenetics.bootstrap_replicates", "replicates", "int", 1000, minimum=0, maximum=100000,
        modes=("all",), category="phylogenetics",
        rationale="Ultrafast bootstrap support; 1000 when computationally feasible.",
        source="Hoang et al. (2018) Mol Biol Evol 35:518-522",
    ),
    # -- resources -----------------------------------------------------------
    Parameter(
        "resources.threads_per_sample", "threads", "int", 2, minimum=1, maximum=1024,
        modes=("all",), category="resources",
        rationale="Conservative default for 8 GB laptops; benchmark before increasing.",
        source="Specification section 15",
    ),
    Parameter(
        "resources.concurrent_samples", "samples", "int", 1, minimum=1, maximum=1000,
        modes=("all",), category="resources",
        rationale="One sample at a time on 8 GB systems.",
        source="Specification section 15",
    ),
    Parameter(
        "resources.mem_mb_per_heavy_rule", "MiB", "int", 6144, minimum=512, maximum=1048576,
        modes=("all",), category="resources",
        rationale="Rule-level memory cap; fail gracefully if peak exceeds capacity.",
        source="Specification section 15",
    ),
    # -- benchmark -----------------------------------------------------------
    Parameter(
        "benchmark.confidence_level", "probability", "float", 0.95, minimum=0.5, maximum=0.999,
        modes=("all",), category="benchmark",
        rationale="Confidence level for Wilson intervals on precision/recall.",
        source="Wilson (1927) JASA 22:209-212",
    ),
    Parameter(
        "benchmark.plan_snp_precision", "fraction", "float", 0.99, minimum=0.0, maximum=1.0,
        modes=("all",), category="benchmark",
        rationale="Planning acceptance target, frozen before external evaluation.",
        source="Specification section 20.9",
    ),
    Parameter(
        "benchmark.plan_snp_recall", "fraction", "float", 0.95, minimum=0.0, maximum=1.0,
        modes=("all",), category="benchmark",
        rationale="Planning acceptance target.",
        source="Specification section 20.9",
    ),
    Parameter(
        "benchmark.plan_indel_precision", "fraction", "float", 0.98, minimum=0.0, maximum=1.0,
        modes=("all",), category="benchmark",
        rationale="Planning acceptance target.",
        source="Specification section 20.9",
    ),
    Parameter(
        "benchmark.plan_indel_recall", "fraction", "float", 0.90, minimum=0.0, maximum=1.0,
        modes=("all",), category="benchmark",
        rationale="Planning acceptance target.",
        source="Specification section 20.9",
    ),
    Parameter(
        "benchmark.plan_lod_sensitivity", "fraction", "float", 0.95, minimum=0.0, maximum=1.0,
        modes=("all",), category="benchmark",
        rationale="Planning LoD sensitivity objective; final LoD may be higher than 5%.",
        source="Specification section 20.9",
    ),
]

REGISTRY: dict[str, Parameter] = {p.name: p for p in _REGISTRY}


def get_parameter(name: str) -> Parameter:
    try:
        return REGISTRY[name]
    except KeyError as exc:  # pragma: no cover - defensive
        raise KeyError(f"Unknown parameter {name!r}; see pfvarx.config.parameters.REGISTRY") from exc


def validate_parameter(name: str, value: Any) -> Any:
    return get_parameter(name).check(value)


def default_for(name: str) -> Any:
    return get_parameter(name).default


def parameters_for_mode(mode: str) -> dict[str, Parameter]:
    return {n: p for n, p in REGISTRY.items() if p.applies_to(mode)}


def registry_records() -> list[dict]:
    return [p.as_dict() for p in _REGISTRY]


def registry_hash() -> str:
    """Stable SHA-256 over the registry contents (order-independent)."""
    import hashlib
    import json

    payload = json.dumps(
        {k: REGISTRY[k].as_dict() for k in sorted(REGISTRY)}, sort_keys=True, default=str
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def unvalidated_parameters() -> list[str]:
    return [n for n, p in REGISTRY.items() if p.policy_status == POLICY_PROVISIONAL]
