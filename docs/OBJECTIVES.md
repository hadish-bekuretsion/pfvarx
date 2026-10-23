# PfVarX — Analysis Objectives

PfVarX is organised around **scientific objectives**, not pipeline stages. You choose the
questions; PfVarX derives the stages, checks the prerequisites, creates only the required
directories, and reports an explicit outcome for every objective you selected.

> **Research use only.** Use only for exploratory research; no scientific or clinical validation is established by this document.

---

## Selecting objectives

```bash
# List them
pixi run -e default pfvarx objectives

# Command line (overrides the configuration file)
pixi run -e default pfvarx run --config cfg.yaml --samplesheet s.tsv --objectives diversity,relatedness

# Or in the configuration
objectives:
  - diversity
  - relatedness
```

Aliases are accepted: `resistance` → `drug_resistance`, `mixed` → `mixed_infection`,
`popstruct` → `population_structure`, `phylo` → `phylogenetics`,
`spatial` → `spatial_temporal`.

Selection order does not matter. PfVarX normalises the set into a documented order so that
two runs with the same objectives produce byte-identical `objective_plan.json` files.

---

## The eight objectives

### `drug_resistance` — Antimalarial drug resistance

**Question:** does this sample carry any variant in a curated resistance marker?

**Requires:** `annotation.gff` (GFF3 gene model) and the bundled resistance catalogue
(17 markers, Pf3D7, catalogue version 1.0.0).

**Artefacts:** `resistance_hits.tsv`, `resistance_summary.json`.

**Scope:** per sample.

**If the GFF is missing** the objective is reported `unavailable` with the missing path
named. It does **not** emit an empty table, because an empty table would be read as "no
resistance mutations found" when in fact nothing was examined.

**Interpretation.** Hits are exploratory annotations, not clinical resistance calls. The
catalogue is marked `NOT_FOR_CLINICAL_USE`.

---

### `mixed_infection` — Mixed infection (minority variants)

**Question:** is there evidence of more than one parasite genotype at a locus?

**Requires:** an evidence table (produced by the caller or imported from VCF).

**Artefacts:** `mixed_infection.json`.

**Scope:** per sample.

**Interpretation.** Reports the distribution of minority-allele bins across callable sites.
A non-zero minority fraction is *evidence consistent with* polyclonality; it is not proof,
because sequencing error and mapping artefacts can mimic it. The caller already models this
in its likelihood; this objective summarises the outcome.

---

### `coi` — Complexity of infection

**Question:** how many parasite clones are present?

**Requires:** an evidence table.

**Artefacts:** `coi.json`.

**Scope:** per sample.

**Interpretation — important.** This objective reports an **ordinal indicator** derived
from the number and fraction of minority variants. It deliberately does **not** report a
numeric COI estimate: the data available to this pipeline cannot support one without
fabricating a number. If you need a numeric COI, use a dedicated method (e.g. a
within-host diversity model fitted to genome-wide allele frequencies) and treat PfVarX's
indicator as corroborating evidence only.

---

### `diversity` — Genetic diversity

**Question:** how diverse is the parasite population within this sample?

**Requires:** an evidence table and the number of callable bases.

**Artefacts:** `diversity.json` (π, segregating-site count, Watterson-type estimate where
computable).

**Scope:** per sample.

**Interpretation.** Diversity is computed over callable sites only. If the callable
fraction is low, π is unreliable regardless of its value; always read it alongside the QC
summary.

---

### `population_structure` — Population structure

**Question:** do samples cluster by population?

**Requires:** a metadata table with a `country` (or `population`) column, and ≥ 2 samples.

**Artefacts:** `population_structure.tsv`, `population_structure.json`.

**Scope:** run.

**Interpretation.** Descriptive aggregation of diversity and allele sharing by group. This
is not a model-based admixture analysis; with small sample sizes, group differences are
likely to be sampling noise.

---

### `relatedness` — Relatedness

**Question:** which samples are related?

**Requires:** ≥ 2 samples with evidence tables.

**Artefacts:** `relatedness.json`, `relatedness.tsv` (pairwise IBS, IBS0, shared and
discordant site counts).

**Scope:** run.

**Interpretation.** Haploid identity-by-state. Pairs sharing an implausibly high fraction of
alleles across many sites are candidates for recent shared ancestry. Relatedness estimates
are only meaningful above a minimum number of shared sites; pairs below that threshold are
reported with their site count so you can judge.

---

### `phylogenetics` — Phylogenetics

**Question:** what is the phylogenetic relationship among samples?

**Requires:** ≥ 2 samples.

**Artefacts:** `tree.nwk` (Newick), `distance_matrix.tsv`.

**Scope:** run.

**Interpretation.** Neighbour-joining tree from the pairwise distance matrix. NJ is fast and
deterministic but assumes no recombination; for recombining *P. falciparum* populations a
network or a recombination-aware model may be more appropriate. The tree is a description of
this dataset, not a species phylogeny.

---

### `spatial_temporal` — Spatial and temporal trends

**Question:** do diversity or allele frequencies vary by place or time?

**Requires:** metadata with `country` **and** `collection_date`, and ≥ 2 samples.

**Artefacts:** `spatial_temporal.tsv`, `temporal_series.tsv`, `spatial_temporal.json`.

**Scope:** run.

**Interpretation.** Descriptive aggregation by country and by collection period. With few
samples per group, apparent trends are not statistically supported; the output includes
group sizes so this is visible.

---

## Prerequisite behaviour

PfVarX checks prerequisites **before** any heavy computation and before creating output
directories. There are two failure modes:

| Situation | Behaviour |
|---|---|
| A prerequisite is missing and the objective cannot run at all | The run **refuses to start**, names the missing item, and suggests the fix. No directories are created. |
| A prerequisite is missing but the run can still proceed | The objective is recorded as `unavailable` with an actionable message; other objectives run normally. |

This distinction matters: `unavailable` means "not examined", not "nothing found".

---

## Dependency expansion

Each objective declares the stages it needs. PfVarX takes the union of the selected
objectives' stages plus the stages implied by each sample's input route, and executes only
those.

Example — `--objectives drug_resistance` on BAM input:

```
route:    alignment_reuse   (no realignment)
stages:   call → annotate
```

Example — `--objectives diversity,mixed_infection,coi` on BAM input:

```
route:    alignment_reuse
stages:   call → diversity → mixed_infection → coi
```

No phylogenetics, no relatedness, no population structure — because none were requested.

---

## Outcomes

Every objective produces exactly one outcome:

| Status | Meaning |
|---|---|
| `completed` | Ran and produced artefacts |
| `failed` | Ran but raised an error; see `error` |
| `skipped` | Deliberately not run |
| `unavailable` | A prerequisite was absent; see `message` |

The run-level status is derived from these: `completed` if all succeeded, `partial_failure`
if any failed, `partial` if any were unavailable, `unavailable` if none could run.

---

## Adding your own objective

1. Add an `Objective` entry to `src/pfvarx/objectives.py` with its `requires`, `artefacts`,
   `needs_metadata`, `needs_catalogue`, `multi_sample` and `scope`.
2. Add a runner to `PER_SAMPLE_RUNNERS` or `RUN_SCOPE_RUNNERS` in
   `src/pfvarx/execution.py`.
3. Add a test in `tests/test_objectives.py`.

The dashboard and exporters pick the new objective up automatically, because they are
generated from the registry and the run manifest rather than from a fixed template.
