# Troubleshooting

This guide covers the errors you are most likely to encounter, what causes them,
and how to resolve them. PfVarX uses a typed exception hierarchy
(`src/pfvarx/errors.py`), so the error class in a traceback tells you which layer
failed.

| Exception | Layer | Typical cause |
|---|---|---|
| `ConfigError` | Configuration | Schema violation, missing required value, precedence conflict |
| `InputValidationError` | Layer 0 | Unreadable, mislabelled, or out-of-scope input |
| `ReferenceMismatchError` | Layer 0/2 | Reference allele does not match the assembly |
| `RoutingError` | Layer 0 | No valid route for a sample |
| `ToolExecutionError` | Layer 1 | External binary failed or is unavailable |
| `ModelError` | Layer 2 | Invalid counts or numerical failure |
| `CalibrationError` | Layer 2 | Calibration data insufficient or malformed |
| `PolicyError` | Layer 2 | Policy status inconsistent with the request |
| `UnsupportedModeError` | Configuration | Mode not in `wgs`/`tbs_no_umi`/`tbs_umi` |
| `OutputError` | Any | Cannot write an output artefact |

---

## 1. Configuration errors

### `ConfigError: extra fields not permitted`

The Pydantic schema uses `extra="forbid"`, so an unknown key is a hard error.
This catches typos that would otherwise be silently ignored.

**Fix.** Compare your config against `config/config.yaml`. Remove or correct the
offending key. The error message names the key and its parent.

### `ConfigError: reference.sha256`

A `VALIDATED` policy requires a real SHA-256. The placeholder
`REQUIRED_SHA256` is rejected.

**Fix.**

```bash
sha256sum reference/Pf3D7.fa
```

Paste the digest into `reference.sha256`. For provisional runs you may leave the
placeholder, but then `ALT_PASS_VALIDATED` is unreachable.

### `ConfigError: reference.assembly_id`

Same as above: a `VALIDATED` policy requires a pinned assembly identifier
(e.g. `Pf3D7_01_v3`).

### `UnsupportedModeError`

The `run.mode` or `--profile` value is not one of `wgs`, `tbs_no_umi`,
`tbs_umi`.

**Fix.** Set `run.mode` to a supported value and use a matching profile.

### Profile merge surprises

Precedence is `defaults < profile < run config < CLI`. If a value seems ignored,
check whether the run config overrides it, or whether the profile sets it after
your expectation. Use `pfvarx validate --json` to see the merged configuration.

---

## 2. Input validation and routing

### `InputValidationError: unsupported or unrecognised file type`

Detection combines suffix and content. A file with a `.bam` suffix that is
actually a text file is rejected.

**Fix.** Confirm the file type with `file <path>` and `samtools quickcheck
<path.bam>`. Rename or re-export as appropriate.

### `InputValidationError: empty file`

A zero-byte file is rejected before any processing.

### `RoutingError: No samples to route`

The samplesheet is empty or has only a header.

**Fix.** Add at least one data row.

### `RoutingError: Duplicate sample IDs after routing`

Two rows resolve to the same `sample_id`.

**Fix.** Make `sample_id` unique.

### `RoutingError: evidence file not found`

A BAM/CRAM path in the samplesheet does not exist relative to the configured
base directory.

**Fix.** Check the path, or pass `--base-dir` to resolve relative paths. Use
`--check-paths` with `validate` to catch this early.

### A reference FASTA is rejected as sample evidence

This is intentional. A reference is not sample evidence and cannot occupy a
sample slot. Supply reads (FASTQ), an alignment (BAM/CRAM), or calls
(VCF/BCF/gVCF).

### A BAM is being realigned

It should not be. If a BAM is realigned, the samplesheet `input_type` is
probably `fastq` or `auto` with a misdetected file. Set `input_type: bam` (or
`cram`) explicitly. The router selects `alignment_reuse` and the alignment rules
are not scheduled.

### A VCF is being re-derived from reads

It should not be. Set `input_type: vcf` (or `bcf`/`gvcf`). The variant-import
route is selected and no calling is performed.

---

## 3. External tool errors

### `ToolExecutionError: bwa-mem2 not found`

The workflow tier is not installed. The Python core does not include the
binaries.

**Fix.** Install and invoke the Pixi `workflow` environment with
`pixi install --locked -e workflow` and `pixi run -e workflow ...`. The legacy
container definitions have not been validated and are not a Pixi-equivalent
deployment. See `docs/INSTALL.md` and `docs/DEPLOYMENT.md`.

### `ToolExecutionError: <tool> exited with code N`

The binary ran but failed. The error includes the command line and stderr tail.

**Fix.** Re-run the command manually to see the full error. Common causes:
missing reference index (`bwa-mem2 index`, `samtools faidx`), disk full, or a
corrupt input.

### `ToolExecutionError: SnpEff database not found`

SnpEff needs a database for the assembly.

**Fix.** Either install/build the database
(`snpEff download <db>` or `snpEff build -gff3`), or set
`annotation.engine: internal` to use the deterministic fallback.

### DEploid / hmmIBD / IQ-TREE wrapper fails closed

These are optional guarded wrappers. They raise when the binary or a required
assumption is absent rather than substituting a weaker method.

**Fix.** Install the tool, or disable the corresponding downstream analysis in
`downstream.enabled`.

---

## 4. Model and calibration errors

### `ModelError: n <= 0 or k < 0 or k > n`

The counts passed to a model are invalid: zero depth, negative count, or more
alternates than depth.

**Fix.** This usually indicates an upstream evidence bug. Check the evidence
table for the site. A zero-depth site should be recorded as `NO_CALL` before
model evaluation, not passed to the model.

### `ModelError: numerical failure`

A likelihood evaluation produced a non-finite value, typically from an extreme
allele fraction or a degenerate beta-binomial parameterisation.

**Fix.** Check `statistics.beta_binomial_rho`; very large `rho` can destabilise
the fit. Use the default `null` (fit from data) or a realistic `≈ 0.001`.

### `CalibrationError`

The control TSV is missing required columns (`depth`, `alt_count`) or has too
few observations to fit.

**Fix.** Provide a TSV with `depth` and `alt_count` columns and enough rows to
estimate the error rate and over-dispersion.

---

## 5. Policy errors

### `ALT_PASS_VALIDATED` never appears

This is by design under `PROVISIONAL_DEVELOPMENT`. The strongest possible call is
`ALT_DETECTED_EXPLORATORY`. To reach `ALT_PASS_VALIDATED` you must:

1. supply an empirical error profile (`statistics.empirical_error_profile`);
2. pin `reference.assembly_id` and `reference.sha256`;
3. set `run.policy_status: VALIDATED`.

See `docs/VALIDATION.md` for what validation requires.

### `PolicyError: posterior threshold required`

Bayesian calling is enabled but `posterior_threshold_if_validated` is unset.

**Fix.** Set an explicit threshold, or disable `statistics.bayesian_enabled`.

---

## 6. Output errors

### `OutputError: cannot write`

The output directory is not writable or the disk is full.

**Fix.** Check permissions and free space. Ensure `outputs.root` is writable.

### Leftover `*.tmp` files

Atomic writes write to a temporary file then rename. A leftover `.tmp` means a
run was killed mid-write. Re-running overwrites it; no data corruption occurs.

### Outputs not overwritten

`outputs.overwrite: false` prevents clobbering existing results. Set it to
`true`, or use a fresh `outputs.root`.

---

## 7. Reporting and dashboard

### PDF export fails

reportlab is not installed.

**Fix.** `pip install -e ".[report]"`.

### Dashboard is blank

The Streamlit app reads a results tree. Point it at a populated `results/`
directory.

**Fix.**

```bash
pfvarx report dashboard --config config/config.yaml --results results
```

### XLSX export fails

openpyxl is missing.

**Fix.** `pip install openpyxl` (included in the core dependencies).

---

## 8. Workflow (Snakemake) issues

### `Samplesheet validation failed`

The Snakefile validates the samplesheet at parse time and aborts on failure. The
error includes the validation result.

**Fix.** Run `pfvarx validate` directly to see the detailed reason.

### Jobs re-run unnecessarily

Snakemake tracks file modification times. If an upstream tool rewrites an
unchanged file, downstream jobs re-run.

**Fix.** Use `--rerun-triggers mtime` or ensure tools do not touch outputs when
nothing changed.

### Out of memory on the local profile

The local profile assumes 8 GB. Large WGS BAMs may need more.

**Fix.** Use the `server` profile or raise `resources.mem_mb`.

### Slurm/Kubernetes jobs fail immediately

The HPC profile contains `REQUIRED_SLURM_ACCOUNT`; the cloud profile assumes a
`pfvarx` namespace and a default storage class.

**Fix.** Edit the profile to match your cluster (`slurm_account`,
`slurm_partition`, `kubernetes-namespace`).

---

## 9. Getting more information

- Run any command with `--json` for machine-readable output.
- Run `pfvarx validate --check-paths --json` to diagnose input problems.
- Run `pfvarx status` to confirm the release status and banner.
- Inspect the evidence table (`<sample>.evidence.tsv`) — it records why depth was
  lost and which filters fired.
- Check the policy identifier (`PVX-...`) to confirm which configuration
  produced a result.
