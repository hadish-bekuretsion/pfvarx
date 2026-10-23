# Parameter Dictionary

PfVarX exposes every analytical parameter through a single typed registry in
`src/pfvarx/config/parameters.py`. The registry is the authoritative source:
the CLI (`pfvarx parameters`), the configuration schema, and this document are
all derived from it. The registry currently contains **57 parameters** across
nine categories, with registry hash `dcf93c05e7a2` (first 12 hex characters).

Each parameter record carries:

- **name** — dotted identifier used in YAML;
- **unit** — physical or logical unit;
- **datatype** — `int`, `float`, `bool`, `str`;
- **default** — proposed default value;
- **minimum / maximum / allowed** — accepted range;
- **modes** — which assay modes the parameter applies to (`wgs`,
  `tbs_no_umi`, `tbs_umi`, or all);
- **rationale** — why the parameter exists;
- **source** — the guideline or publication the default is drawn from;
- **policy_status** — `PROVISIONAL_DEVELOPMENT`, `VALIDATED`, or `DERIVED`;
- **calibration_dataset_id** — dataset used to calibrate (null when provisional);
- **version** — registry version.

Every parameter in the shipped registry is `PROVISIONAL_DEVELOPMENT`. This is a
deliberate honesty constraint: none of the defaults has been calibrated against
an empirical truth set in this release, so none may be presented as validated.
A parameter becomes `VALIDATED` only after calibration against a documented
control dataset, at which point `calibration_dataset_id` is populated.

Registry presence and schema acceptance do **not** guarantee that a field changes
every execution path. Some values are planning metadata or are not yet wired to
the active caller/workflow. The README's [field-by-field configuration reference](../README.md#field-by-field-reference)
marks the runtime status of each accepted field; use it before changing a
parameter. In particular, do not treat `require_checksums`, `fastp_max_n_fraction`,
`empirical_error_profile`, benchmark target values, or output retention/atomic
write preferences as active runtime controls in this release.

---

## 1. Using the registry

```bash
# Print the full registry as JSON
pixi run -e default pfvarx parameters

# Filter to one mode
pixi run -e default pfvarx parameters --mode wgs

# Export to a file
pixi run -e default pfvarx parameters --out parameters.json
```

Run Python snippets inside the Pixi `default` environment. Programmatically:

```python
from pfvarx.config.parameters import REGISTRY, get_parameter, registry_hash

p = get_parameter("statistics.beta_binomial_rho")
print(p.default, p.minimum, p.maximum, p.policy_status)
print(registry_hash()[:12])
```

Validation of a proposed value:

```python
from pfvarx.config.parameters import validate_parameter
validate_parameter("statistics.exploratory_fdr_q", 0.05)   # ok
validate_parameter("statistics.exploratory_fdr_q", 1.5)    # raises
```

---

## 2. Scientific justification

The defaults are engineering starting points anchored to established practice,
not tuned cutoffs. The principal sources are:

- **MalariaGEN Pf7 bioinformatics methods** — the community reference pipeline
  for *P. falciparum* WGS, which informs mapping quality, base quality, and
  depth conventions.
- **Gardner et al. (2002), Nature 419:498-511** — the Pf3D7 reference genome
  and the haploid nuclear genome assumption (`reference.nuclear_ploidy = 1`).
- **Chen et al. (2018), Bioinformatics 34:884-890** — the fastp tool and its
  quality/length filtering conventions.
- **Li (2013), arXiv:1303.3997** — BWA-MEM alignment and MAPQ semantics.
- `exploratory_fdr_q` and `sensitivity_fdr_q` are legacy names: the caller
  applies them as decision cutoffs to Bonferroni-adjusted p-values, not as
  Benjamini–Hochberg FDR levels. The `bh_fdr` helper is not the caller's path.
- **Self & Liang (1987)** — boundary chi-square mixture for the LRT null.
- **SAMtools/CRAM specification and GA4GH refget** — reference checksum
  conventions (`reference.sha256`).

The full citation list with verified bibliographic details is in
`docs/METHODS.md`.

---

## 3. Parameter tables

The tables below are generated directly from the registry. The **Range** column
shows the accepted interval or the allowed enumeration; **Modes** shows `all`
when the parameter applies to every assay; **Status** is the policy status.


### Reference (3 parameters)

| Parameter | Unit | Type | Default | Range | Modes | Status | Source |
|---|---|---|---|---|---|---|---|
| `reference.assembly_id` | identifier | str | `REQUIRED_PINNED_ASSEMBLY_ID` |  | all | PROVISIONAL_DEVELOPMENT | MalariaGEN Pf7 bioinformatics methods (2023) |
| `reference.sha256` | hex digest | str | `REQUIRED_SHA256` |  | all | PROVISIONAL_DEVELOPMENT | SAMtools/CRAM specification; GA4GH refget |
| `reference.nuclear_ploidy` | copies | int | `1` | 1–1 | all | PROVISIONAL_DEVELOPMENT | Gardner et al. (2002) Nature 419:498-511 |

### Preprocessing (8 parameters)

| Parameter | Unit | Type | Default | Range | Modes | Status | Source |
|---|---|---|---|---|---|---|---|
| `preprocessing.fastp_min_phred` | Phred | int | `20` | 0–40 | all | PROVISIONAL_DEVELOPMENT | fastp documentation (Chen et al. 2018, Bioinformatics 34:884-890) |
| `preprocessing.fastp_min_length` | bp | int | `50` | 0–1000 | all | PROVISIONAL_DEVELOPMENT | MalariaGEN Pf7 methods; BWA-MEM2 recommendations |
| `preprocessing.fastp_max_n_fraction` | fraction | float | `0.05` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | fastp documentation |
| `preprocessing.mapping_quality_min` | MAPQ | int | `30` | 0–60 | all | PROVISIONAL_DEVELOPMENT | Li (2013) arXiv:1303.3997; MalariaGEN Pf7 methods |
| `preprocessing.base_quality_min` | Phred | int | `20` | 0–40 | all | PROVISIONAL_DEVELOPMENT | MalariaGEN Pf7 methods |
| `preprocessing.count_overlapping_mates_once` | boolean | bool | `True` |  | all | PROVISIONAL_DEVELOPMENT | samtools mpileup -x / GATK read-overlap conventions |
| `preprocessing.wgs_mark_duplicates` | boolean | bool | `True` |  | wgs | PROVISIONAL_DEVELOPMENT | Picard MarkDuplicates documentation |
| `preprocessing.duplicate_exclusion_from_evidence` | boolean | bool | `True` |  | wgs | PROVISIONAL_DEVELOPMENT | MalariaGEN Pf7 methods |

### Callability (7 parameters)

| Parameter | Unit | Type | Default | Range | Modes | Status | Source |
|---|---|---|---|---|---|---|---|
| `callability.wgs_min_independent_depth` | fragments | int | `10` | 1–100000 | wgs | PROVISIONAL_DEVELOPMENT | MalariaGEN Pf7 methods; Pf7 QC description (2023) |
| `callability.targeted_min_independent_depth` | fragments | int | `100` | 1–1000000 | tbs_no_umi | PROVISIONAL_DEVELOPMENT | Reporting considerations for targeted malaria molecular surveillance (2025) |
| `callability.umi_min_consensus_molecules_majority` | molecules | int | `20` | 1–100000 | tbs_umi | PROVISIONAL_DEVELOPMENT | fgbio documentation; UMI consensus literature |
| `callability.umi_min_consensus_molecules_exploratory` | molecules | int | `100` | 1–1000000 | tbs_umi | PROVISIONAL_DEVELOPMENT | fgbio documentation |
| `callability.wgs_high_depth_multiplier_for_review` | multiplier | float | `5.0` | 1.0–100.0 | wgs | PROVISIONAL_DEVELOPMENT | MalariaGEN Pf7 methods |
| `callability.wgs_max_depth_floor` | fragments | int | `500` | 0–100000 | wgs | PROVISIONAL_DEVELOPMENT | MalariaGEN Pf7 methods |
| `callability.wgs_min_callable_fraction` | fraction | float | `0.8` | 0.0–1.0 | wgs | PROVISIONAL_DEVELOPMENT | MalariaGEN Pf7 methods |

### Candidate (6 parameters)

| Parameter | Unit | Type | Default | Range | Modes | Status | Source |
|---|---|---|---|---|---|---|---|
| `candidate.snp_min_alt_observations` | fragments | int | `2` | 1–100 | all | PROVISIONAL_DEVELOPMENT | Development policy; cf. BCFtools mpileup default |
| `candidate.indel_min_alt_observations` | fragments | int | `3` | 1–100 | all | PROVISIONAL_DEVELOPMENT | Development policy; cf. FreeBayes/GATK indel defaults |
| `candidate.exploratory_min_alt_fraction` | fraction | float | `0.01` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Development policy |
| `candidate.max_small_indel_bp` | bp | int | `50` | 1–50 | all | PROVISIONAL_DEVELOPMENT | Specification section 2 |
| `candidate.homopolymer_review_length` | bp | int | `6` | 1–50 | all | PROVISIONAL_DEVELOPMENT | Development policy; homopolymer indel literature |
| `candidate.strand_bias_alpha` | probability | float | `0.01` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Guo et al. (2012) strand-bias assessment |

### Statistics (12 parameters)

| Parameter | Unit | Type | Default | Range | Modes | Status | Source |
|---|---|---|---|---|---|---|---|
| `statistics.count_model` | enum | str | `auto` | one of: binomial, beta_binomial, auto | all | PROVISIONAL_DEVELOPMENT | Development policy; beta-binomial overdispersion literature |
| `statistics.error_rate` | probability | float | `None` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Development policy; requires assay-specific control data |
| `statistics.beta_binomial_rho` | fraction | float | `None` | 0.0–0.999999 | all | PROVISIONAL_DEVELOPMENT | Development policy |
| `statistics.min_alt_fraction_h1` | fraction | float | `0.01` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Development policy; must be replaced by validated LoD |
| `statistics.lrt_null_calibration` | enum | str | `bootstrap` | one of: bootstrap, chi2_mixture, none | all | PROVISIONAL_DEVELOPMENT | Self & Liang (1987); parametric bootstrap literature |
| `statistics.bootstrap_replicates` | replicates | int | `200` | 0–100000 | all | PROVISIONAL_DEVELOPMENT | Development policy |
| `statistics.exploratory_fdr_q` | decision cutoff | float | `0.05` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Development policy; applied after Bonferroni correction |
| `statistics.sensitivity_fdr_q` | decision cutoff | float | `0.01` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Development policy; applied after Bonferroni correction |
| `statistics.bayesian_enabled` | boolean | bool | `False` |  | all | PROVISIONAL_DEVELOPMENT | Specification section 6.2 Model D |
| `statistics.bayesian_prior_h1` | probability | float | `0.01` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Development policy |
| `statistics.posterior_threshold_if_validated` | probability | float | `None` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Development policy |
| `statistics.low_depth_fallback` | enum | str | `conservative` | one of: conservative, fail | all | PROVISIONAL_DEVELOPMENT | Specification section 20.5 |

### Reporting (5 parameters)

| Parameter | Unit | Type | Default | Range | Modes | Status | Source |
|---|---|---|---|---|---|---|---|
| `reporting.minority_candidate_fraction` | fraction | float | `0.05` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Development policy |
| `reporting.intermediate_fraction` | fraction | float | `0.2` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Development policy |
| `reporting.predominant_fraction` | fraction | float | `0.8` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Development policy |
| `reporting.pass_requires_assay_validated_lod` | boolean | bool | `True` |  | all | PROVISIONAL_DEVELOPMENT | Specification section 20.5 |
| `reporting.uncertainty_always_visible` | boolean | bool | `True` |  | all | PROVISIONAL_DEVELOPMENT | Specification section 20.8 |

### Population (6 parameters)

| Parameter | Unit | Type | Default | Range | Modes | Status | Source |
|---|---|---|---|---|---|---|---|
| `population.sample_missingness_max` | fraction | float | `0.2` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Development policy |
| `population.variant_site_missingness_max` | fraction | float | `0.1` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Development policy |
| `population.maf_min` | fraction | float | `0.05` | 0.0–0.5 | all | PROVISIONAL_DEVELOPMENT | Development policy |
| `population.ld_window_snps` | SNPs | int | `50` | 1–100000 | all | PROVISIONAL_DEVELOPMENT | PLINK --indep-pairwise convention |
| `population.ld_step_snps` | SNPs | int | `5` | 1–10000 | all | PROVISIONAL_DEVELOPMENT | PLINK --indep-pairwise convention |
| `population.ld_r2_max` | r^2 | float | `0.2` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | PLINK --indep-pairwise convention |

### Phylogenetics (1 parameters)

| Parameter | Unit | Type | Default | Range | Modes | Status | Source |
|---|---|---|---|---|---|---|---|
| `phylogenetics.bootstrap_replicates` | replicates | int | `1000` | 0–100000 | all | PROVISIONAL_DEVELOPMENT | Hoang et al. (2018) Mol Biol Evol 35:518-522 |

### Resources (3 parameters)

| Parameter | Unit | Type | Default | Range | Modes | Status | Source |
|---|---|---|---|---|---|---|---|
| `resources.threads_per_sample` | threads | int | `2` | 1–1024 | all | PROVISIONAL_DEVELOPMENT | Specification section 15 |
| `resources.concurrent_samples` | samples | int | `1` | 1–1000 | all | PROVISIONAL_DEVELOPMENT | Specification section 15 |
| `resources.mem_mb_per_heavy_rule` | MiB | int | `6144` | 512–1048576 | all | PROVISIONAL_DEVELOPMENT | Specification section 15 |

### Benchmark (6 parameters)

| Parameter | Unit | Type | Default | Range | Modes | Status | Source |
|---|---|---|---|---|---|---|---|
| `benchmark.confidence_level` | probability | float | `0.95` | 0.5–0.999 | all | PROVISIONAL_DEVELOPMENT | Wilson (1927) JASA 22:209-212 |
| `benchmark.plan_snp_precision` | fraction | float | `0.99` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Specification section 20.9 |
| `benchmark.plan_snp_recall` | fraction | float | `0.95` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Specification section 20.9 |
| `benchmark.plan_indel_precision` | fraction | float | `0.98` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Specification section 20.9 |
| `benchmark.plan_indel_recall` | fraction | float | `0.9` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Specification section 20.9 |
| `benchmark.plan_lod_sensitivity` | fraction | float | `0.95` | 0.0–1.0 | all | PROVISIONAL_DEVELOPMENT | Specification section 20.9 |

---

## 4. Notes on selected parameters

**`statistics.error_rate`** and **`statistics.beta_binomial_rho`** are treated by
the caller as configured fixed null parameters. The caller does not represent
uncertainty from their estimation or context-specific error rates. `rho → 0`
recovers the binomial model. Calibration utilities can estimate parameters, but
this does not propagate their estimation uncertainty into per-call inference or
establish an assay LoD.

**`statistics.lrt_null_calibration`** and **`statistics.bootstrap_replicates`**
are retained for compatibility with standalone LRT helper routines. They do not
change the current caller's final `lrt_pvalue`, which is the exact one-sided
count-tail p-value. The caller records a note when bootstrap was requested; it
does not use that bootstrap result for its decision path.

**`statistics.exploratory_fdr_q`** and **`sensitivity_fdr_q`** retain legacy
`fdr_q` names, but are cutoffs (`0.05` and `0.01`) applied to the Bonferroni
family-adjusted count-tail p-value by the caller. They do not configure BH/FDR.
The caller corrects across the complete pre-filter reference-scan upper bound,
including all possible configured anchored insertions/deletions up to
`candidate.max_small_indel_bp`; this can be extremely conservative and reduce
power. This correction is per sample, not a cohort-wide adjustment across
multiple samples. `bh_fdr` remains a helper, not the caller's inference path.

**`callability.wgs_min_independent_depth`** is the minimum informative depth
for a WGS site to be callable. Below it, a hard gate forces `NO_CALL` with
`LowEffectiveDepth`. The targeted default (`100`) reflects the higher depth of
targeted sequencing.

**`candidate.max_small_indel_bp`** caps INDEL size at 50 bp, enforcing the
specification's small-INDEL scope. Larger events are rejected at evidence
extraction.

**`reporting.*_fraction`** define the minority bins: `0.05` (minority
candidate), `0.20` (intermediate), `0.80` (predominant). These are reporting
boundaries, not calling thresholds.

**`resources.*`** control workflow parallelism and are consumed by the
Snakemake profiles, not by the statistical caller.

---

## 5. Changing parameters safely

1. Edit the run configuration (or a profile), not the registry source.
2. Validate the change: `pfvarx validate --config <cfg> --samplesheet <tsv>`.
3. Re-run and compare the evidence table against the previous run.
4. Record the change and its justification in your run notes; the
   `policy_id` in every output encodes the configuration hash, so the change is
   traceable.

Do not present any result produced under `PROVISIONAL_DEVELOPMENT` parameters as
validated. See `docs/VALIDATION.md` for what validation would require.
