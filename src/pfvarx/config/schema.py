"""Pydantic v2 configuration schema for PfVarX.

The schema enforces structural validity, numeric bounds, mutually exclusive
sequencing modes and conditional requirements (for example ``tbs_umi`` demands
a verified UMI design and consensus settings). It deliberately does *not*
claim scientific validity: ``policy_status`` records whether the numeric
values are provisional or validated.
"""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from pfvarx.config import parameters as P
from pfvarx.errors import ConfigError

Mode = Literal["wgs", "tbs_no_umi", "tbs_umi"]
PolicyStatus = Literal["PROVISIONAL_DEVELOPMENT", "VALIDATED"]
Intent = Literal["research", "surveillance", "clinical_research"]
InputType = Literal["auto", "fastq", "bam", "cram", "vcf", "bcf", "gvcf"]

_REQUIRED_PLACEHOLDER = re.compile(r"^REQUIRED_")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


class StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True, str_strip_whitespace=True)


class RunConfig(StrictModel):
    mode: Mode
    intent: Intent = "research"
    policy_status: PolicyStatus = "PROVISIONAL_DEVELOPMENT"
    run_id: str | None = None
    description: str = ""

    @field_validator("run_id")
    @classmethod
    def _safe_run_id(cls, v: str | None) -> str | None:
        if v is None:
            return v
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}", v):
            raise ValueError(
                "run_id must be 1-64 characters of letters, digits, dot, underscore or hyphen "
                "and must start with an alphanumeric character"
            )
        return v


class ReferenceConfig(StrictModel):
    assembly_id: str
    fasta: str
    sha256: str
    callable_mask_bed: str | None = None
    snpeff_db: str | None = None
    snpeff_jar: str | None = None

    @field_validator("sha256")
    @classmethod
    def _sha(cls, v: str) -> str:
        if _REQUIRED_PLACEHOLDER.match(v):
            return v
        if not _SHA256.match(v.lower()):
            raise ValueError("reference.sha256 must be a 64-character hexadecimal SHA-256 digest")
        return v.lower()


class InputsConfig(StrictModel):
    samplesheet: str
    root: str = "data"
    raw_dir: str | None = None
    prepared_dir: str | None = None
    metadata_dir: str | None = None
    input_type: InputType = "auto"
    accept_existing_bam: bool = True
    require_checksums: bool = False


class PreprocessingConfig(StrictModel):
    fastp_min_phred: int = Field(20, ge=0, le=40)
    fastp_min_length: int = Field(50, ge=0, le=1000)
    fastp_max_n_fraction: float = Field(0.05, ge=0.0, le=1.0)
    adapter_autodetect: bool = True
    primer_bed: str | None = None
    mapping_quality_min: int = Field(30, ge=0, le=60)
    base_quality_min: int = Field(20, ge=0, le=40)
    count_overlapping_mates_once: bool = True
    wgs_mark_duplicates: bool = True
    duplicate_exclusion_from_evidence: bool = True
    umi_design: str | None = None
    umi_min_family_size: int = Field(2, ge=1, le=100)


class CallabilityConfig(StrictModel):
    wgs_min_independent_depth: int = Field(10, ge=1)
    targeted_min_independent_depth: int = Field(100, ge=1)
    umi_min_consensus_molecules_majority: int = Field(20, ge=1)
    umi_min_consensus_molecules_exploratory: int = Field(100, ge=1)
    wgs_high_depth_multiplier_for_review: float = Field(5.0, ge=1.0)
    wgs_max_depth_floor: int = Field(500, ge=0)
    wgs_min_callable_fraction: float = Field(0.80, ge=0.0, le=1.0)

    def min_depth_for_mode(self, mode: str) -> int:
        if mode == "wgs":
            return self.wgs_min_independent_depth
        if mode == "tbs_no_umi":
            return self.targeted_min_independent_depth
        return self.umi_min_consensus_molecules_majority


class CandidateConfig(StrictModel):
    snp_min_alt_observations: int = Field(2, ge=1)
    indel_min_alt_observations: int = Field(3, ge=1)
    exploratory_min_alt_fraction: float = Field(0.01, ge=0.0, le=1.0)
    max_small_indel_bp: int = Field(50, ge=1, le=50)
    homopolymer_review_length: int = Field(6, ge=1, le=50)
    strand_bias_alpha: float = Field(0.01, ge=0.0, le=1.0)


class StatisticsConfig(StrictModel):
    count_model: Literal["binomial", "beta_binomial", "auto"] = "auto"
    empirical_error_profile: str | None = None
    error_rate: float | None = Field(None, ge=0.0, le=1.0)
    beta_binomial_rho: float | None = Field(None, ge=0.0, lt=1.0)
    min_alt_fraction_h1: float = Field(0.01, ge=0.0, le=1.0)
    lrt_null_calibration: Literal["bootstrap", "chi2_mixture", "none"] = "bootstrap"
    bootstrap_replicates: int = Field(200, ge=0, le=100000)
    exploratory_fdr_q: float = Field(0.05, ge=0.0, le=1.0)
    sensitivity_fdr_q: float = Field(0.01, ge=0.0, le=1.0)
    bayesian_enabled: bool = False
    bayesian_prior_h1: float = Field(0.01, ge=0.0, le=1.0)
    posterior_threshold_if_validated: float | None = Field(None, ge=0.0, le=1.0)
    low_depth_fallback: Literal["conservative", "fail"] = "conservative"

    @model_validator(mode="after")
    def _rho_and_bayes(self) -> StatisticsConfig:
        if self.beta_binomial_rho is not None and not (0.0 <= self.beta_binomial_rho < 1.0):
            raise ValueError("beta_binomial_rho must satisfy 0 <= rho < 1")
        if self.bayesian_enabled and self.posterior_threshold_if_validated is None:
            raise ValueError(
                "bayesian_enabled=true requires posterior_threshold_if_validated "
                "(a calibrated gate must be stated explicitly)"
            )
        return self


class ReportingConfig(StrictModel):
    minority_candidate_fraction: float = Field(0.05, ge=0.0, le=1.0)
    intermediate_fraction: float = Field(0.20, ge=0.0, le=1.0)
    predominant_fraction: float = Field(0.80, ge=0.0, le=1.0)
    pass_requires_assay_validated_lod: bool = True
    uncertainty_always_visible: bool = True

    @model_validator(mode="after")
    def _ordered(self) -> ReportingConfig:
        if not (
            self.minority_candidate_fraction
            <= self.intermediate_fraction
            <= self.predominant_fraction
        ):
            raise ValueError("reporting fraction bins must be non-decreasing")
        return self


class BenchmarkConfig(StrictModel):
    confidence_level: float = Field(0.95, ge=0.5, le=0.999)
    plan_snp_precision: float = Field(0.99, ge=0.0, le=1.0)
    plan_snp_recall: float = Field(0.95, ge=0.0, le=1.0)
    plan_indel_precision: float = Field(0.98, ge=0.0, le=1.0)
    plan_indel_recall: float = Field(0.90, ge=0.0, le=1.0)
    plan_lod_sensitivity: float = Field(0.95, ge=0.0, le=1.0)


class ResourcesConfig(StrictModel):
    threads_per_sample: int = Field(2, ge=1, le=1024)
    concurrent_samples: int = Field(1, ge=1, le=1000)
    mem_mb_per_heavy_rule: int = Field(6144, ge=512, le=1048576)


class OutputsConfig(StrictModel):
    root: str = "results"
    logs_dir: str | None = None
    work_dir: str | None = None
    retain_raw_evidence: bool = True
    atomic_writes: bool = True
    overwrite: bool = False


class AnnotationConfig(StrictModel):
    enabled: bool = True
    engine: Literal["internal", "auto"] = "auto"
    snpeff_db: str | None = None
    #: GFF3 gene model used to derive amino-acid changes.  Without it the
    #: drug-resistance objective reports ``unavailable`` rather than an empty
    #: (and therefore misleading) resistance table.
    gff: str | None = None


class DownstreamConfig(StrictModel):
    enabled: list[
        Literal[
            "drug_resistance",
            "mixed_infection",
            "coi",
            "diversity",
            "population_structure",
            "relatedness",
            "phylogenetics",
            "spatial_temporal",
        ]
    ] = Field(default_factory=list)
    resistance_catalogue: str | None = None
    metadata_table: str | None = None


class PfVarXConfig(StrictModel):
    schema_version: int = Field(1, ge=1, le=1)
    project: str = "PfVarX"
    run: RunConfig
    reference: ReferenceConfig
    inputs: InputsConfig
    preprocessing: PreprocessingConfig = Field(default_factory=PreprocessingConfig)
    callability: CallabilityConfig = Field(default_factory=CallabilityConfig)
    candidate: CandidateConfig = Field(default_factory=CandidateConfig)
    statistics: StatisticsConfig = Field(default_factory=StatisticsConfig)
    reporting: ReportingConfig = Field(default_factory=ReportingConfig)
    benchmark: BenchmarkConfig = Field(default_factory=BenchmarkConfig)
    resources: ResourcesConfig = Field(default_factory=ResourcesConfig)
    outputs: OutputsConfig = Field(default_factory=OutputsConfig)
    annotation: AnnotationConfig = Field(default_factory=AnnotationConfig)
    downstream: DownstreamConfig = Field(default_factory=DownstreamConfig)
    #: Analysis objectives the user explicitly asked for.  When empty, the
    #: legacy ``downstream.enabled`` list is used instead so that older
    #: configuration files keep working unchanged.
    objectives: list[str] = Field(default_factory=list)

    # -- objectives ----------------------------------------------------------
    @field_validator("objectives")
    @classmethod
    def _known_objectives(cls, value: list[str]) -> list[str]:
        from pfvarx.objectives import normalise

        resolved: list[str] = []
        for name in value:
            canonical = normalise(name)
            if canonical not in resolved:
                resolved.append(canonical)
        return resolved

    @model_validator(mode="after")
    def _objectives_mirror_downstream(self) -> PfVarXConfig:
        """Keep ``objectives`` and ``downstream.enabled`` consistent.

        ``objectives`` is the user-facing surface; ``downstream.enabled`` is the
        internal list consumed by routing.  Whichever was supplied populates the
        other, so there is exactly one effective selection.
        """
        if self.objectives and not self.downstream.enabled:
            self.downstream.enabled = list(self.objectives)
        elif self.downstream.enabled and not self.objectives:
            self.objectives = list(self.downstream.enabled)
        elif (
            self.objectives
            and self.downstream.enabled
            and list(self.objectives) != list(self.downstream.enabled)
        ):
            raise ConfigError(
                "objectives and downstream.enabled disagree",
                cause=f"objectives={self.objectives} downstream.enabled={self.downstream.enabled}",
                action="Set only `objectives:` at the top level and remove downstream.enabled",
            )
        return self

    # -- conditional validity ------------------------------------------------
    @model_validator(mode="after")
    def _mode_conditions(self) -> PfVarXConfig:
        if self.run.mode == "tbs_umi":
            if not self.preprocessing.umi_design:
                raise ConfigError(
                    "tbs_umi requires preprocessing.umi_design describing a verified UMI layout",
                    cause="UMI design not supplied",
                    action="Provide a validated UMI design identifier or choose tbs_no_umi",
                    consequence="UMI mode is refused rather than silently downgraded",
                )
            if self.preprocessing.umi_min_family_size < 2:
                raise ConfigError("tbs_umi requires umi_min_family_size >= 2")
        if self.run.policy_status == "VALIDATED":
            if _REQUIRED_PLACEHOLDER.match(self.reference.assembly_id):
                raise ConfigError("VALIDATED policy requires a pinned reference.assembly_id")
            if self.statistics.error_rate is None:
                raise ConfigError(
                    "VALIDATED policy requires statistics.error_rate from an empirical profile",
                    action="Supply a signed calibration profile",
                )
            if self.run.mode == "tbs_umi" and self.statistics.beta_binomial_rho is None:
                raise ConfigError("VALIDATED tbs_umi policy requires a calibrated beta_binomial_rho")
        if self.reporting.pass_requires_assay_validated_lod and self.run.policy_status != "VALIDATED":
            # Not an error: it forces the policy engine to emit RESEARCH_CANDIDATE.
            pass
        if self.run.intent == "clinical_research" and self.run.policy_status != "VALIDATED":
            raise ConfigError(
                "clinical_research intent requires policy_status=VALIDATED",
                action="Use intent=research for provisional development runs",
            )
        return self

    # -- helpers -------------------------------------------------------------
    def is_placeholder(self, field: str) -> bool:
        value = self.reference.assembly_id if field == "assembly_id" else getattr(self.reference, field, "")
        return bool(value and _REQUIRED_PLACEHOLDER.match(str(value)))

    def unresolved_placeholders(self) -> list[str]:
        out = []
        if _REQUIRED_PLACEHOLDER.match(self.reference.assembly_id):
            out.append("reference.assembly_id")
        if _REQUIRED_PLACEHOLDER.match(self.reference.sha256):
            out.append("reference.sha256")
        if _REQUIRED_PLACEHOLDER.match(self.statistics.empirical_error_profile or ""):
            out.append("statistics.empirical_error_profile")
        return out

    def canonical_dict(self, exclude_run_id: bool = True) -> dict:
        d = self.model_dump(mode="json")
        if exclude_run_id:
            d.get("run", {}).pop("run_id", None)
        return d

    def config_hash(self, exclude_run_id: bool = True) -> str:
        payload = json.dumps(self.canonical_dict(exclude_run_id), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode()).hexdigest()

    def policy_id(self) -> str:
        """Immutable identifier combining config hash and parameter-registry hash."""
        return f"PVX-{self.run.mode}-{self.run.policy_status[:4]}-{self.config_hash()[:12]}-{P.registry_hash()[:8]}"

    def effective_parameters(self) -> dict:
        """Flatten the effective scientific settings with provenance."""
        flat = {
            "preprocessing.fastp_min_phred": self.preprocessing.fastp_min_phred,
            "preprocessing.fastp_min_length": self.preprocessing.fastp_min_length,
            "preprocessing.fastp_max_n_fraction": self.preprocessing.fastp_max_n_fraction,
            "preprocessing.mapping_quality_min": self.preprocessing.mapping_quality_min,
            "preprocessing.base_quality_min": self.preprocessing.base_quality_min,
            "preprocessing.count_overlapping_mates_once": self.preprocessing.count_overlapping_mates_once,
            "preprocessing.wgs_mark_duplicates": self.preprocessing.wgs_mark_duplicates,
            "preprocessing.duplicate_exclusion_from_evidence": self.preprocessing.duplicate_exclusion_from_evidence,
            "callability.wgs_min_independent_depth": self.callability.wgs_min_independent_depth,
            "callability.targeted_min_independent_depth": self.callability.targeted_min_independent_depth,
            "callability.umi_min_consensus_molecules_majority": self.callability.umi_min_consensus_molecules_majority,
            "callability.umi_min_consensus_molecules_exploratory": self.callability.umi_min_consensus_molecules_exploratory,
            "callability.wgs_high_depth_multiplier_for_review": self.callability.wgs_high_depth_multiplier_for_review,
            "callability.wgs_max_depth_floor": self.callability.wgs_max_depth_floor,
            "candidate.snp_min_alt_observations": self.candidate.snp_min_alt_observations,
            "candidate.indel_min_alt_observations": self.candidate.indel_min_alt_observations,
            "candidate.exploratory_min_alt_fraction": self.candidate.exploratory_min_alt_fraction,
            "candidate.max_small_indel_bp": self.candidate.max_small_indel_bp,
            "candidate.homopolymer_review_length": self.candidate.homopolymer_review_length,
            "candidate.strand_bias_alpha": self.candidate.strand_bias_alpha,
            "statistics.count_model": self.statistics.count_model,
            "statistics.error_rate": self.statistics.error_rate,
            "statistics.beta_binomial_rho": self.statistics.beta_binomial_rho,
            "statistics.min_alt_fraction_h1": self.statistics.min_alt_fraction_h1,
            "statistics.lrt_null_calibration": self.statistics.lrt_null_calibration,
            "statistics.bootstrap_replicates": self.statistics.bootstrap_replicates,
            "statistics.exploratory_fdr_q": self.statistics.exploratory_fdr_q,
            "statistics.sensitivity_fdr_q": self.statistics.sensitivity_fdr_q,
            "statistics.bayesian_enabled": self.statistics.bayesian_enabled,
            "statistics.bayesian_prior_h1": self.statistics.bayesian_prior_h1,
            "reporting.minority_candidate_fraction": self.reporting.minority_candidate_fraction,
            "reporting.intermediate_fraction": self.reporting.intermediate_fraction,
            "reporting.predominant_fraction": self.reporting.predominant_fraction,
            "resources.threads_per_sample": self.resources.threads_per_sample,
            "resources.concurrent_samples": self.resources.concurrent_samples,
            "resources.mem_mb_per_heavy_rule": self.resources.mem_mb_per_heavy_rule,
        }
        out = {}
        for name, value in flat.items():
            meta = P.REGISTRY.get(name)
            out[name] = {
                "value": value,
                "unit": meta.unit if meta else None,
                "policy_status": meta.policy_status if meta else None,
                "source": meta.source if meta else None,
                "version": meta.version if meta else None,
                "applicable": meta.applies_to(self.run.mode) if meta else None,
            }
        return out

    def scientific_warnings(self) -> list[str]:
        warns: list[str] = []
        if self.run.policy_status == "PROVISIONAL_DEVELOPMENT":
            warns.append(
                "policy_status=PROVISIONAL_DEVELOPMENT: no call may be labelled ALT_PASS_VALIDATED; "
                "outputs are RESEARCH_CANDIDATE/UNVALIDATED."
            )
        if self.statistics.error_rate is None:
            warns.append(
                "statistics.error_rate is unset; the caller will run sensitivity scenarios "
                "p0 in {0.001, 0.005, 0.01} and mark results UNVALIDATED."
            )
        if self.statistics.beta_binomial_rho is None:
            warns.append("beta_binomial_rho is unset; overdispersion is not calibrated (OVERDISPERSION_NOT_CALIBRATED).")
        if self.run.mode == "tbs_no_umi":
            warns.append(
                "tbs_no_umi: PCR/amplicon duplication is not equivalent to WGS duplication; "
                "duplicate marking is disabled and assay-specific calibration is required."
            )
        if self.run.intent != "research":
            warns.append("intent is not 'research'; confirm institutional governance before use.")
        return warns

    def validate_paths(self, base_dir: Path | None = None, require_reference: bool = True) -> list[str]:
        """Check declared files exist. Returns warnings for optional files."""
        base = Path(base_dir) if base_dir else Path.cwd()
        warnings: list[str] = []
        required = [("reference.fasta", self.reference.fasta), ("inputs.samplesheet", self.inputs.samplesheet)]
        for label, rel in required:
            p = Path(rel)
            p = p if p.is_absolute() else base / p
            if not p.exists():
                raise ConfigError(
                    f"Required file not found for {label}: {p}",
                    cause="Missing input file",
                    action=f"Create or correct {label}",
                    consequence="Run cannot start",
                )
        for label, rel in [
            ("reference.callable_mask_bed", self.reference.callable_mask_bed),
            ("preprocessing.primer_bed", self.preprocessing.primer_bed),
            ("downstream.resistance_catalogue", self.downstream.resistance_catalogue),
            ("downstream.metadata_table", self.downstream.metadata_table),
        ]:
            if rel:
                p = Path(rel)
                p = p if p.is_absolute() else base / p
                if not p.exists():
                    warnings.append(f"{label} not found: {p}")
        return warnings
