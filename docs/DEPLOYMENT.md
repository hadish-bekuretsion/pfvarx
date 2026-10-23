# Deployment

The repository contains local and shared-server profiles plus Slurm and
Kubernetes profile templates. Only the local profile has been smoke-tested here,
using a synthetic paired-FASTQ dataset; no cluster or container deployment has
been validated. This document describes the templates, resource assumptions,
and operational practices for reproducible research runs.

---

## 1. Deployment tiers

| Tier | Executor | Typical use | Minimum resources |
|---|---|---|---|
| Local | `local` | Development, single sample, demo | 8 GB RAM, 2–4 cores |
| Server | `local` | Small cohorts on a shared machine | 64 GB RAM, 32 cores |
| HPC | `slurm` profile template | Slurm executor plugin is installed in the workflow environment; cluster execution is not validated | Cluster-specific |
| Cloud | `kubernetes` profile template | Kubernetes executor plugin is installed in the workflow environment; cluster execution is not validated | Cluster-specific |

The Python core (caller, QC, reporting, benchmark) runs on any of these without
Snakemake. Snakemake is required only for multi-stage, multi-sample
orchestration.

For a Pixi-managed local install, prepare the workflow environment once and
prefix the profile commands below with `pixi run -e workflow`:

```bash
pixi install --locked -e workflow
pixi run -e workflow snakemake --profile workflow/profiles/local
```

The `workflow` environment contains Snakemake, workflow binaries, test tools,
and the PyPI packages `snakemake-executor-plugin-slurm` 2.8.x and
`snakemake-executor-plugin-kubernetes` 0.5.x. Their executor interfaces are
compatible with the locked Snakemake 9.x environment. `pixi.lock` records
resolved builds for `linux-64`, `linux-aarch64`, and `osx-64`; use
`pixi install --locked` to verify that the manifest and lock agree. The plugins
are installed locally, but no Slurm or Kubernetes cluster was available for
validation, and package resolution does not establish scientific
reproducibility.

---

## 2. Local workstation

```bash
pixi run -e workflow snakemake --profile workflow/profiles/local
```

Profile (`workflow/profiles/local/config.yaml`):

```yaml
executor: local
jobs: 1
cores: 4
resources:
  mem_mb: 8192
  runtime: 720
latency-wait: 30
rerun-incomplete: true
keep-going: true
```

Suitable for a laptop or a small VM. The 8 GB figure is the platform's stated
minimum; whole-genome BAMs from deep sequencing may need more, in which case use
the server profile.

---

## 3. Shared server

```bash
pixi run -e workflow snakemake --profile workflow/profiles/server
```

Profile:

```yaml
executor: local
jobs: 8
cores: 32
resources:
  mem_mb: 65536
  runtime: 1440
latency-wait: 60
rerun-incomplete: true
keep-going: true
```

Runs up to eight jobs concurrently on a single machine. Use this when you have
a dedicated analysis server with ample RAM.

---

## 4. HPC (Slurm) — profile template; cluster validation still required

The repository includes `workflow/profiles/hpc/config.yaml`. The workflow
environment now includes the PyPI distribution
[`snakemake-executor-plugin-slurm`](https://pypi.org/project/snakemake-executor-plugin-slurm/)
2.8.x, which requires Snakemake executor interfaces 9.4.x and is locked with
Snakemake 9.x. The profile is still a template, not a verified Slurm
deployment. Configure the account and partition, provide Slurm client commands
and cluster access, then test it on the target cluster before use.

Profile:

```yaml
executor: slurm
jobs: 100
cores: 64
default-resources:
  slurm_partition: cpu
  slurm_account: REQUIRED_SLURM_ACCOUNT
  mem_mb: 16384
  runtime: 2880
latency-wait: 120
rerun-incomplete: true
keep-going: true
```

The template also contains cluster-specific account and partition placeholders.
Editing those values alone does not provide cluster credentials, scheduler
access, or prove the profile operational.

For large cohorts, stage inputs on a shared filesystem visible to all nodes, and
write outputs to a scratch filesystem, then copy the final results tree to
long-term storage.

---

## 5. Cloud (Kubernetes) — profile template; cluster validation still required

The repository includes `workflow/profiles/cloud/config.yaml`. The workflow
environment includes the PyPI distribution
[`snakemake-executor-plugin-kubernetes`](https://pypi.org/project/snakemake-executor-plugin-kubernetes/)
0.5.x, whose executor interface requirements are compatible with the locked
Snakemake 9.x environment. This profile is not verified as an operational
deployment. Validate cluster credentials, namespace, pod/container execution,
and input/output storage against a real cluster before use.

Profile:

```yaml
executor: kubernetes
jobs: 50
cores: 32
default-resources:
  mem_mb: 16384
  runtime: 2880
kubernetes-namespace: pfvarx
latency-wait: 120
rerun-incomplete: true
keep-going: true
```

The profile's namespace and storage settings need cluster-specific validation.
For example:

```bash
kubectl create namespace pfvarx
kubectl get storageclass
```

Mount inputs and outputs via a persistent volume claim. The profile assumes the
`pfvarx` namespace; change it if your cluster uses a different convention.

---

## 6. Containers

**Status: unverified legacy recipes; not part of the verified Pixi deployment path.**
The current Dockerfile and Apptainer definition still build from
`environment.yml` with micromamba rather than installing from `pixi.lock`.
They have not been built or tested in this audit and must not be treated as a
reproducible production deployment. Use the Pixi environments on a supported
host. The container definitions require migration to Pixi and independent build
and runtime tests before this status can change.

### 6.1 Docker

```bash
docker build -t pfvarx:0.9.0 -f container/Dockerfile .
docker run --rm -v "$PWD:/data" pfvarx:0.9.0 status
docker run --rm -v "$PWD:/data" pfvarx:0.9.0 \
    demo --config /opt/pfvarx/config/profiles/wgs.yaml --out /data/demo
```

These commands are retained for reference only; they describe the legacy
micromamba-based image and do not satisfy the current Pixi-lock release gate.

### 6.2 Apptainer / Singularity

```bash
apptainer build pfvarx_0.9.0.sif container/pfvarx.def
apptainer run pfvarx_0.9.0.sif status
apptainer exec --bind "$PWD:/data" pfvarx_0.9.0.sif \
    pfvarx demo --config /opt/pfvarx/config/profiles/wgs.yaml --out /data/demo
```

The definition file includes a `%test` section, but it has not been built or
executed. Do not treat its presence as container validation or use it as a
recommended HPC route until it is migrated to Pixi and tested.

---

## 7. Resource planning

Resource requirements scale with cohort size and sequencing depth, not with the
number of variants. The dominant costs are:

- **Alignment** — memory scales with reference index size and thread count;
  BWA-MEM2 uses roughly 1–2 GB per thread for the *P. falciparum* genome.
- **Duplicate marking** — Picard holds a sort buffer proportional to the number
  of reads; 16 GB is comfortable for typical WGS BAMs.
- **Calling** — the pure-Python caller streams the BAM region by region; memory
  is modest (hundreds of MB) but throughput depends on candidate density.
- **Annotation** — the internal GFF3 consequence implementation is lightweight;
  SnpEff is installed but not integrated into per-call annotation.

`resources.threads_per_sample` and `resources.mem_mb_per_heavy_rule` are used
where the workflow rules request them. The accepted `concurrent_samples` value
does not override the shipped profile's fixed `jobs` setting. These resource
settings are separate from statistical parameters and do not affect call
accuracy.

---

## 8. Storage layout

Recommended production layout:

```
/data/
├── reference/          # pinned FASTA + indexes (read-only)
├── inputs/             # FASTQ/BAM/VCF (read-only)
├── results/            # outputs (read-write)
├── work/               # Snakemake intermediate files (scratch)
└── logs/               # per-job logs
```

Keep `reference/` and `inputs/` immutable and checksummed. Point the configured
`outputs.work_dir` at fast scratch storage and copy final results to long-term
storage. Verify every profile-specific storage and path setting on the target
execution environment before use.

---

## 9. Reproducibility in production

1. **Pin the reference.** Record `reference.assembly_id` and
   `reference.sha256`. Never run a `VALIDATED` policy against an unpinned
   assembly.
2. **Pin the environment.** Keep `pixi.toml` and `pixi.lock` together, use
   `pixi install --locked`, and review/commit lockfile changes with dependency
   changes. The lock pins software builds for declared platforms; it does not
   pin external reference data or assay calibration.
3. **Pin the configuration.** Commit the run config. The `policy_id` in every
   output encodes the configuration hash, so a result can always be traced to
   its config.
4. **Seed randomised steps.** Pass `--seed` to `benchmark` and `demo`; the
   bootstrap takes a seed from the configuration.
5. **Keep the evidence table.** It is the scientific record; the VCF alone is
   not sufficient for audit.

---

## 10. Monitoring and logging

- Snakemake writes per-job logs; the profiles set `printshellcmds: true` for the
  local profile.
- Every run writes `results/run_qc.json`, which aggregates QC metrics and the
  overall status.
- `pfvarx status` reports the release-status matrix at any time.
- For long runs, use `--keep-going` (set in all profiles) so independent samples
  continue after a single failure, and inspect the failed job's log.

---

## 11. Security and data governance

PfVarX is a research tool. It does not implement authentication, encryption at
rest, or access control; those are the responsibility of the hosting
environment. For human-genomics-adjacent data (e.g. linked clinical metadata),
apply your institution's data-governance policy. The platform never transmits
data externally; all processing is local or within your cluster.

---

## 12. Upgrade procedure

1. Pull the new version and update the selected Pixi environment
   (`pixi install --locked -e workflow`).
2. Re-run `pixi run -e default status` and `pixi run -e default test` to confirm the install.
3. Re-run `pixi run -e default pfvarx validate` on a representative samplesheet.
4. Compare a small test sample's evidence table against the previous version.
   Because outputs are idempotent, any difference is attributable to a code or
   parameter change and should be reviewed.
5. Update the container image tag.

---

## 13. Deployment checklist

- [ ] Reference indexed and checksummed
- [ ] Samplesheet validated (`pixi run -e default pfvarx validate --check-paths`)
- [ ] Configuration reviewed; placeholders replaced
- [ ] Local/server profile selected; HPC/Kubernetes profiles used only after cluster-specific validation
- [ ] Environment resolution recorded (`pixi.lock` or image tag)
- [ ] Output and scratch directories writable
- [ ] `pixi run -e default pfvarx status` shows the expected banner
- [ ] Test run on one sample succeeds end to end
- [ ] QC summary reviewed
- [ ] Evidence tables retained
