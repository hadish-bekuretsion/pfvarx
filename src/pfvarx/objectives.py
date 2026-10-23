"""Analysis objectives and their dependency graph.

An *objective* is a scientific question the user wants answered, not a stage of
the pipeline.  This module is the single source of truth for

* which objectives exist and what each one means,
* which stages each objective requires,
* which objectives need auxiliary inputs (metadata table, resistance
  catalogue), and
* which artefacts each objective produces.

The workflow is built from the selected objectives, so a run that asks only for
drug resistance never touches the diversity, relatedness or phylogenetics code
and never creates directories for them.  This is the mechanism that satisfies
the requirement that the pipeline "must not automatically execute every
available analysis".
"""

from __future__ import annotations

from dataclasses import dataclass, field

from pfvarx.errors import ConfigError

# -- objective identifiers --------------------------------------------------
OBJ_DRUG_RESISTANCE = "drug_resistance"
OBJ_MIXED_INFECTION = "mixed_infection"
OBJ_COI = "coi"
OBJ_DIVERSITY = "diversity"
OBJ_POPULATION_STRUCTURE = "population_structure"
OBJ_RELATEDNESS = "relatedness"
OBJ_PHYLOGENETICS = "phylogenetics"
OBJ_SPATIAL_TEMPORAL = "spatial_temporal"

ALL_OBJECTIVES = (
    OBJ_DRUG_RESISTANCE,
    OBJ_MIXED_INFECTION,
    OBJ_COI,
    OBJ_DIVERSITY,
    OBJ_POPULATION_STRUCTURE,
    OBJ_RELATEDNESS,
    OBJ_PHYLOGENETICS,
    OBJ_SPATIAL_TEMPORAL,
)

# Objectives that describe more than one sample.  Selecting any of these makes
# the run multi-sample: the planner will not silently compute a pairwise metric
# from a single sample.
MULTI_SAMPLE_OBJECTIVES = frozenset(
    {OBJ_POPULATION_STRUCTURE, OBJ_RELATEDNESS, OBJ_PHYLOGENETICS}
)

# Objectives that read a sample metadata table (one row per sample).
METADATA_OBJECTIVES = frozenset({OBJ_POPULATION_STRUCTURE, OBJ_SPATIAL_TEMPORAL})


@dataclass(frozen=True)
class Objective:
    """A scientific objective and the machinery it needs."""

    id: str
    label: str
    summary: str
    #: Stages that must run to satisfy the objective.  ``call`` implies
    #: candidate discovery and evidence extraction, which are never optional
    #: when variants are needed.
    requires: tuple[str, ...]
    #: Files written into ``<root>/objectives/<id>/``.
    artefacts: tuple[str, ...]
    #: Whether a per-sample metadata row is required.
    needs_metadata: bool = False
    #: Whether a resistance catalogue is required.
    needs_catalogue: bool = False
    #: Whether the analysis is meaningful with a single sample.
    multi_sample: bool = False
    #: Whether the result is a per-sample table or a run-level table.
    scope: str = "per_sample"

    def as_dict(self) -> dict:
        return {
            "id": self.id,
            "label": self.label,
            "summary": self.summary,
            "requires": list(self.requires),
            "artefacts": list(self.artefacts),
            "needs_metadata": self.needs_metadata,
            "needs_catalogue": self.needs_catalogue,
            "multi_sample": self.multi_sample,
            "scope": self.scope,
        }


OBJECTIVES: dict[str, Objective] = {
    OBJ_DRUG_RESISTANCE: Objective(
        id=OBJ_DRUG_RESISTANCE,
        label="Antimalarial drug resistance",
        summary=(
            "Annotate called variants against a curated resistance-marker "
            "catalogue and report, per sample, which known or candidate "
            "resistance substitutions are supported."
        ),
        requires=("call", "annotate"),
        artefacts=("resistance_hits.tsv", "resistance_summary.json"),
        needs_catalogue=True,
        scope="per_sample",
    ),
    OBJ_MIXED_INFECTION: Objective(
        id=OBJ_MIXED_INFECTION,
        label="Mixed infection (minority variants)",
        summary=(
            "Summarise within-host variant fractions and classify sites into "
            "minority, intermediate and predominant bins."
        ),
        requires=("call",),
        artefacts=("mixed_infection.json",),
        scope="per_sample",
    ),
    OBJ_COI: Objective(
        id=OBJ_COI,
        label="Complexity of infection (COI)",
        summary=(
            "Report the complexity-of-infection estimate when a validated "
            "estimator is available, and state explicitly when it is not."
        ),
        requires=("call",),
        artefacts=("coi.json",),
        scope="per_sample",
    ),
    OBJ_DIVERSITY: Objective(
        id=OBJ_DIVERSITY,
        label="Genetic diversity",
        summary=(
            "Per-sample nucleotide diversity, expected heterozygosity and the "
            "count of polymorphic sites."
        ),
        requires=("call",),
        artefacts=("diversity.json",),
        scope="per_sample",
    ),
    OBJ_POPULATION_STRUCTURE: Objective(
        id=OBJ_POPULATION_STRUCTURE,
        label="Population structure",
        summary=(
            "Group samples by metadata fields and summarise diversity and "
            "variant burden per group."
        ),
        requires=("call",),
        artefacts=("population_structure.tsv", "population_structure.json"),
        needs_metadata=True,
        scope="run",
    ),
    OBJ_RELATEDNESS: Objective(
        id=OBJ_RELATEDNESS,
        label="Relatedness",
        summary=(
            "Pairwise identity-by-state between samples from shared, "
            "sufficiently covered sites."
        ),
        requires=("call",),
        artefacts=("relatedness.json", "relatedness.tsv"),
        multi_sample=True,
        scope="run",
    ),
    OBJ_PHYLOGENETICS: Objective(
        id=OBJ_PHYLOGENETICS,
        label="Phylogenetics",
        summary=(
            "Neighbour-joining tree from the pairwise distance matrix, "
            "written as Newick."
        ),
        requires=("call",),
        artefacts=("tree.nwk", "distance_matrix.tsv"),
        multi_sample=True,
        scope="run",
    ),
    OBJ_SPATIAL_TEMPORAL: Objective(
        id=OBJ_SPATIAL_TEMPORAL,
        label="Spatial and temporal trends",
        summary=(
            "Aggregate per-sample summaries by location and collection date "
            "using the metadata table."
        ),
        requires=("call",),
        artefacts=("spatial_temporal.tsv", "spatial_temporal.json"),
        needs_metadata=True,
        scope="run",
    ),
}

#: Human-readable aliases accepted on the command line and in configuration.
ALIASES = {
    "resistance": OBJ_DRUG_RESISTANCE,
    "drug-resistance": OBJ_DRUG_RESISTANCE,
    "drugresistance": OBJ_DRUG_RESISTANCE,
    "mixed": OBJ_MIXED_INFECTION,
    "mixed-infection": OBJ_MIXED_INFECTION,
    "minority": OBJ_MIXED_INFECTION,
    "complexity": OBJ_COI,
    "complexity_of_infection": OBJ_COI,
    "population-structure": OBJ_POPULATION_STRUCTURE,
    "popstruct": OBJ_POPULATION_STRUCTURE,
    "structure": OBJ_POPULATION_STRUCTURE,
    "related": OBJ_RELATEDNESS,
    "ibd": OBJ_RELATEDNESS,
    "phylo": OBJ_PHYLOGENETICS,
    "tree": OBJ_PHYLOGENETICS,
    "spatial": OBJ_SPATIAL_TEMPORAL,
    "temporal": OBJ_SPATIAL_TEMPORAL,
    "spatiotemporal": OBJ_SPATIAL_TEMPORAL,
}


@dataclass
class ResolvedObjectives:
    """The validated, dependency-expanded objective set for a run."""

    selected: list[str] = field(default_factory=list)
    stages: list[str] = field(default_factory=list)
    needs_metadata: bool = False
    needs_catalogue: bool = False
    multi_sample: bool = False
    warnings: list[str] = field(default_factory=list)

    @property
    def is_empty(self) -> bool:
        return not self.selected

    def as_dict(self) -> dict:
        return {
            "selected": list(self.selected),
            "labels": [OBJECTIVES[o].label for o in self.selected],
            "stages": list(self.stages),
            "needs_metadata": self.needs_metadata,
            "needs_catalogue": self.needs_catalogue,
            "multi_sample": self.multi_sample,
            "warnings": list(self.warnings),
        }


def normalise(name: str) -> str:
    """Map a user-supplied objective name to its canonical identifier."""
    key = str(name).strip().lower().replace(" ", "_")
    if key in OBJECTIVES:
        return key
    if key in ALIASES:
        return ALIASES[key]
    raise ConfigError(
        f"Unknown analysis objective: {name!r}",
        cause=f"Known objectives: {', '.join(ALL_OBJECTIVES)}",
        action=(
            "Use one of the identifiers above, or run 'pfvarx objectives' to "
            "list them with descriptions."
        ),
    )


def resolve(
    selected: list[str] | tuple[str, ...] | None,
    *,
    n_samples: int = 1,
    have_metadata: bool = False,
    have_catalogue: bool = False,
) -> ResolvedObjectives:
    """Validate the selection and expand it to the required stage set.

    Raises :class:`ConfigError` when the selection is unknown or when an
    objective needs an input that was not supplied.  Failing here, before any
    computation, is deliberate: a run must not start and then discover halfway
    through that the metadata table it needs does not exist.
    """
    if not selected:
        return ResolvedObjectives()

    canonical: list[str] = []
    for name in selected:
        oid = normalise(name)
        if oid not in canonical:
            canonical.append(oid)

    # Preserve the documented order rather than the user's typing order, so
    # that two runs with the same set produce byte-identical manifests.
    ordered = [o for o in ALL_OBJECTIVES if o in canonical]

    warnings: list[str] = []
    stages: list[str] = []
    needs_metadata = False
    needs_catalogue = False
    multi_sample = False

    for oid in ordered:
        obj = OBJECTIVES[oid]
        for stage in obj.requires:
            if stage not in stages:
                stages.append(stage)
        if obj.needs_metadata:
            needs_metadata = True
        if obj.needs_catalogue:
            needs_catalogue = True
        if obj.multi_sample:
            multi_sample = True

    if needs_metadata and not have_metadata:
        missing = [OBJECTIVES[o].id for o in ordered if OBJECTIVES[o].needs_metadata]
        raise ConfigError(
            "The selected objectives require a sample metadata table, but none "
            f"was supplied: {', '.join(missing)}",
            cause="downstream.metadata_table is unset and no --metadata was given",
            action=(
                "Provide a TSV/CSV with one row per sample and a 'sample_id' "
                "column, then set downstream.metadata_table or pass --metadata."
            ),
        )

    if needs_catalogue and not have_catalogue:
        raise ConfigError(
            "The drug-resistance objective requires a resistance catalogue, but "
            "none could be loaded",
            cause="annotation catalogue path is unset and package data is missing",
            action=(
                "Set downstream.resistance_catalogue to a YAML catalogue, or "
                "reinstall the package so its bundled catalogue is present."
            ),
        )

    if multi_sample and n_samples < 2:
        offenders = [OBJECTIVES[o].label for o in ordered if OBJECTIVES[o].multi_sample]
        raise ConfigError(
            "These objectives compare samples and need at least two samples: "
            f"{', '.join(offenders)}",
            cause=f"Only {n_samples} sample(s) were supplied",
            action=(
                "Add more samples to the samplesheet, or deselect the "
                "multi-sample objectives."
            ),
        )

    # A single-sample run that also asks for a run-level objective is fine as
    # long as that objective is not multi-sample; note it rather than fail.
    if n_samples == 1 and any(OBJECTIVES[o].scope == "run" for o in ordered):
        warnings.append(
            "A run-level objective was selected with a single sample; its "
            "output will be a one-row table."
        )

    return ResolvedObjectives(
        selected=ordered,
        stages=stages,
        needs_metadata=needs_metadata,
        needs_catalogue=needs_catalogue,
        multi_sample=multi_sample,
        warnings=warnings,
    )


def describe() -> list[dict]:
    """Machine-readable catalogue of objectives, for ``pfvarx objectives``."""
    return [OBJECTIVES[o].as_dict() for o in ALL_OBJECTIVES]


def format_table() -> str:
    """Human-readable objective listing for the CLI."""
    lines = [
        f"{'ID':<22} {'LABEL':<38} NEEDS",
        "-" * 92,
    ]
    for oid in ALL_OBJECTIVES:
        obj = OBJECTIVES[oid]
        needs = []
        if obj.needs_metadata:
            needs.append("metadata")
        if obj.needs_catalogue:
            needs.append("catalogue")
        if obj.multi_sample:
            needs.append(">=2 samples")
        lines.append(
            f"{obj.id:<22} {obj.label:<38} {', '.join(needs) if needs else '-'}"
        )
    lines.append("")
    lines.append("Select with:  --objectives drug_resistance,diversity")
    lines.append("Or in config: objectives: [drug_resistance, diversity]")
    return "\n".join(lines)
