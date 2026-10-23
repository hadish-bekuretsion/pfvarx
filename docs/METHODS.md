# Analytical Methods

This document describes the statistical and computational methods implemented in
PfVarX, with the rationale for each choice and citations to the primary
literature. It is written to be read alongside `docs/PARAMETERS.md` (which lists
every parameter) and `docs/VALIDATION.md` (which states what has and has not
been validated).

---

## 1. Scope and assumptions

PfVarX targets **single-nucleotide polymorphisms and small insertions/deletions
(1–50 bp)** in *Plasmodium falciparum*. It assumes:

- a **haploid** nuclear genome (`reference.nuclear_ploidy = 1`), consistent with
  the reference genome assembly (Gardner et al. 2002);
- **Illumina paired-end short-read** data, either whole-genome (WGS) or
  targeted (TBS), with or without unique molecular identifiers (UMIs);
- the possibility of **polyclonal infections**, in which a single sample may
  contain multiple parasite genotypes, so that the alternate allele fraction at
  a site can lie anywhere in `(0, 1)` rather than being restricted to `0` or `1`;
- a **pinned reference assembly** whose identity and checksum are recorded.

Structural variants, copy-number variation, large INDELs (>50 bp), and
non-Illumina platforms are explicitly out of scope. The caller contains **no
machine-learning component**; every decision is a transparent statistical rule.

---

## 2. Input validation and routing

Input files are typed by suffix and content. Content sniffing checks BGZF magic
bytes for compressed formats and the `##fileformat` header for VCF/BCF, so a
mislabelled file is rejected rather than misprocessed. Checksums are computed
with SHA-256 (National Institute of Standards and Technology 2015) and recorded
in a provenance record alongside each artefact.

Routing honours two hard constraints from the specification:

1. **A valid BAM/CRAM is never realigned.** If a sample supplies an aligned
   file, the router selects `alignment_reuse` and the alignment rules are not
   scheduled. Re-alignment is both wasteful and a source of irreproducibility.
2. **A VCF/BCF/gVCF follows the variant-import route.** Pre-called variants are
   imported with `model = "imported"` and `calibration_profile = "IMPORTED"`.
   They are never re-derived from reads, and they can never reach
   `ALT_PASS_VALIDATED`, because the underlying evidence was not observed by
   PfVarX.

A reference FASTA is recognised as *not sample evidence* and cannot occupy a
sample slot.

---

## 3. Read preprocessing

### 3.1 Trimming and quality filtering

Adapter sequences and low-quality tails are removed with fastp (Chen et al.
2018). Defaults are a minimum Phred of 20, a minimum retained length of 50 bp,
and a maximum ambiguous-base fraction of 0.05. These are conventional Illumina
thresholds; the exact values are provisional and configurable.

### 3.2 Alignment

Reads are aligned with BWA-MEM2 (Li 2013; Vasimuddin et al. 2019), which
implements the BWA-MEM seed-and-extend algorithm with a faster, SIMD-accelerated
backend. Mapping quality (MAPQ) follows the BWA-MEM definition; reads below
`preprocessing.mapping_quality_min` (default 30) are excluded from evidence.
The MAPQ 30 convention is the standard "uniquely mapped" cutoff and is used by
the MalariaGEN Pf7 pipeline (MalariaGEN 2023).

### 3.3 Duplicate handling

For WGS, PCR duplicates are marked with Picard MarkDuplicates (Broad Institute).
Duplicate fragments are excluded from evidence by default
(`preprocessing.duplicate_exclusion_from_evidence`). For targeted sequencing
without UMIs, duplicate marking is disabled by default because the fragment
distribution is not PCR-dominated in the same way; this is a provisional choice
and is documented as such.

### 3.4 UMI consensus

For TBS with UMIs, read families are collapsed with fgbio
(`fgbio GroupReadsByUmi` / `CallMolecularConsensusReads`). A UMI design
identifier is **required** (`umi_design: REQUIRED_UMI_DESIGN_ID`); the pipeline
fails closed if it is absent. The minimum family size defaults to 2. UMI
consensus reduces both PCR and sequencing errors by requiring independent
molecules to agree (Kivioja et al. 2011; Fu et al. 2011; Schmitt et al. 2012;
Kinde et al. 2011; Newman et al. 2016).

---

## 4. Evidence extraction

Evidence extraction is a CIGAR-aware pileup implemented in
`calling/evidence.py`. For each candidate site it records:

- **reference, alternate, and other counts** after filtering;
- **informative depth** — the number of reads that actually inform the
  reference-versus-alternate decision;
- **allele fraction** `alt_count / informative_depth` with a **Jeffreys
  credible interval** (Jeffreys 1946);
- **mean base quality** for reference and alternate observations;
- **strand counts** — forward and reverse observations of each allele;
- **mean read position** and a **read-position bias** statistic;
- a **Fisher exact test** p-value for strand bias (Fisher 1935);
- **homopolymer run length** at the site;
- **left-alignment status** for INDELs.

Filtering rules:

- overlapping mate pairs are counted once
  (`preprocessing.count_overlapping_mates_once`), preventing double-counting of
  a single DNA fragment;
- low-MAPQ reads, low-base-quality bases, secondary alignments, and duplicate
  fragments are excluded and counted separately, so the evidence table records
  *why* depth was lost;
- INDELs are represented in anchored form and left-aligned to a canonical
  position, so the same event is not counted twice under different
  representations.

The Jeffreys interval is used rather than a Wald interval because the latter
behaves poorly near the boundaries `f = 0` and `f = 1`, which are exactly the
regions of interest for minority-variant detection (Brown, Cai & DasGupta 2001).

---

## 5. Statistical models

PfVarX implements four count models. All likelihoods are evaluated in log space
using `scipy.special.gammaln`, `betaln`, and `logsumexp` to avoid underflow at
high depth.

### 5.1 Model A — binomial

Under the null hypothesis `H0`, alternate observations arise from sequencing
error at rate `e`:

```
H0:  k ~ Binomial(n, e)
H1:  k ~ Binomial(n, p1(f)),   f >= configured AF boundary
p1(f) = f + (1-f)*e
```

The log-likelihood is

```
logL = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1) + k*log(p) + (n-k)*log(1-p)
```

The binomial model assumes independent, identically distributed observations.
It is the simplest model and understates variance when over-dispersion is
present.

### 5.2 Model B — beta-binomial

To accommodate over-dispersion (extra-binomial variance from library
preparation, mapping ambiguity, and within-sample heterogeneity), the
alternate-observation probability is drawn from a Beta distribution:

```
p ~ Beta(alpha, beta)
alpha = mu * (1 - rho) / rho
beta  = (1 - mu) * (1 - rho) / rho
```

where `mu` is the mean error (or allele) fraction and `rho` is the
over-dispersion parameter. As `rho → 0`, the beta-binomial converges to the
binomial. The log-likelihood uses the beta-binomial coefficient:

```
logL = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
     + betaln(k+alpha, n-k+beta) - betaln(alpha, beta)
```

The beta-binomial is an option for modelling over-dispersion; it is not a
universally validated default. Its `rho` must be estimated from appropriate
assay-specific controls before inferential use. The code treats the configured
`rho` as fixed and does not propagate its estimation uncertainty. A literature
value or the synthetic benchmark setting is not an assay calibration.

### 5.3 Model C — likelihoods and caller inference

`count_models.evaluate_hypotheses` compares the error-only null with a **composite
alternative**: true allele fraction (AF) at or above the configured boundary
`statistics.min_alt_fraction_h1`. The alternative likelihood is maximized over
the allowed AF interval; it is not a point alternative. The reported
`lrt_statistic = 2 * (logL1 - logL0)` describes that likelihood comparison.

This LRT is **not the caller's inferential p-value path**. For each candidate,
`calling.policy` uses the exact one-sided upper tail of the ALT count under the
configured binomial or beta-binomial null. That p-value is written to the legacy
field `lrt_pvalue`; its name does not mean that it is an LRT-calibrated p-value.
Bootstrap and boundary-mixture LRT routines remain in `count_models` as helpers
(the analytic boundary-mixture approximation follows Self & Liang 1987), but
they do not replace this caller decision rule. `bh_fdr` likewise remains a
model helper and is not the caller's multiple-testing procedure.

The caller multiplies each raw count-tail p-value by a **Bonferroni family-size
upper bound** and caps the result at one. The bound is calculated over the
complete requested reference scan before evidence/candidate filtering and
`max_candidates` truncation. Per scanned anchor it counts three SNP alternatives
plus every possible anchored insertion and deletion sequence up to
`candidate.max_small_indel_bp`. This is an intentionally broad bound and can be
extremely conservative, especially with the configured small-indel universe;
reduced power is a material limitation. The configured `exploratory_fdr_q` and
`sensitivity_fdr_q` values are used as decision cutoffs after this family-wise
adjustment despite their legacy `fdr_q` names. The evidence field
`lrt_pvalue_adjusted` holds the Bonferroni-adjusted value, not a BH q-value.
The family adjustment is performed per sample; it does not provide a cohort-wide
family-wise guarantee across multiple samples.

The tail probability is exact **conditional on the configured null parameters**
and informative depth. `statistics.error_rate` and `statistics.beta_binomial_rho`
are treated as fixed values; uncertainty from their estimation and context-
specific error variation are not represented by this caller test. An absent
error rate triggers a provisional placeholder for exploratory reporting, not an
empirical calibration.

### 5.4 Model D — exploratory Bayesian quantity

The optional Bayesian calculation is disabled by default and requires an explicit
`posterior_threshold_if_validated`. When enabled it uses the configured AF boundary and prior without integrating
over an allele-fraction prior; it is a **point-boundary-only exploratory**
quantity, not an assay-calibrated posterior for biological
truth, and is never sufficient to validate or certify a call. Uncertainty in the
configured error rate or `rho` is not propagated.

### 5.5 What is *not* done

A raw LRT statistic is **never** relabelled as a posterior probability, and an
uncalibrated metric is **never** encoded as a VCF `QUAL`. The writer requires a
`VALIDATED` policy status, a posterior, and `ALT_PASS_VALIDATED`; because the
current caller cannot produce a validated PASS without assay-validation/LoD
provenance, the normal caller path leaves `QUAL` missing. This prevents a
likelihood ratio from being presented as a Phred-scaled confidence.

---

## 6. Decision policy

The decision policy maps evidence and model output to one of five categories:

| Category | Meaning |
|---|---|
| `NO_CALL` | A hard gate failed; no decision is possible |
| `REFERENCE_SUPPORTED` | Evidence is consistent with the reference allele |
| `ALT_DETECTED_EXPLORATORY` | Alternate is supported but not validated |
| `ALT_PASS_VALIDATED` | Reserved category; the current caller cannot emit it without verifiable assay-validation/LoD provenance |
| `ALT_UNCERTAIN` | Evidence exists but is flagged for review |

**Hard gates** force `NO_CALL`:

- `LowEffectiveDepth` — informative depth below the callability threshold;
- `OutsideCallableRegion` — the site is outside the callable mask;
- `ReferenceMismatch` — the reference allele does not match the assembly.

**Review filters** downgrade to `ALT_UNCERTAIN`:

- `StrandBiasReview` — Fisher exact strand-bias p-value below
  `candidate.strand_bias_alpha`;
- `HomopolymerIndelReview` — INDEL in a homopolymer run at or above
  `candidate.homopolymer_review_length`;
- `AmbiguousMapping`, `UMIModelInvalid`, `UncalibratedModel` — contextual flags.

The policy status and configuration schema alone do not establish assay
validation. The current caller consumes no verifiable assay-validation/LoD
provenance object and explicitly downgrades any would-be validated PASS;
therefore it cannot produce a validated `ALT_PASS_VALIDATED` call. Calls are at
most exploratory (or uncertain/no-call according to evidence and filters), even
if a run is labelled `VALIDATED`. Tests and synthetic demos establish software
behaviour only, not biological validation.

**Minority bins** classify the allele fraction for reporting: `f < 0.01` is
below threshold; `0.01–0.05` is `LOW_FRACTION_EXPLORATORY`; `0.05–0.20` is
`MINORITY_CANDIDATE`; `0.20–0.80` is `INTERMEDIATE_FRACTION`; `≥ 0.80` is
`PREDOMINANT_ALT`. The intermediate bin is the signal used for mixed-infection
inference.

---

## 7. Output schema

Each site yields a record in a **46-column evidence table** and a corresponding
VCF record. The evidence table carries the reference/alternate/uninformative
counts, depth, allele fraction and interval, model identity and version, both
log-likelihoods, the LRT statistic and p-values, the exploratory posterior quantity (when
present), the calibration-profile label, the filter status and reasons, uncertainty
flags, sequencing mode, and UMI status. The VCF carries the decision category
and policy identifier as `PVX_DECISION` and `PVX_POLICY_ID`.

This separation is deliberate. The VCF is the interchange format; the evidence
table is the scientific record. A consumer who needs only calls reads the VCF; a
consumer who needs to audit or re-analyse reads the evidence table.

---

## 8. Annotation

Per-call annotation uses the internal GFF3 + codon-translation engine only.
`annotation.engine` accepts `internal` or `auto`; selecting SnpEff is rejected by
the configuration schema, and the per-call path does not invoke SnpEff. The
internal engine classifies consequences when the needed GFF3 and FASTA inputs are
available; otherwise consequences can be reported as unavailable/unannotated.
Resistance markers are matched against a curated catalogue
(`config/catalogue/pf_resistance.yaml`); matches are evidence annotations, not
clinical interpretations. Unit tests establish implementation behaviour, not
scientific accuracy of annotations.

---

## 9. Downstream analyses

### 9.1 Diversity

Nucleotide diversity `π` is the average pairwise difference per site (Nei & Li
1979). Watterson's theta `θ_W` is estimated from the number of segregating sites
(Watterson 1975). Tajima's `D` compares the two (Tajima 1989). All are computed
over callable bases, with missingness thresholds from the parameter registry.

### 9.2 Mixed infection

Samples are classified by the number of variants in the intermediate
allele-fraction bin (`0.20–0.80`). Three or more intermediate variants yields
`MIXED_INFECTION_LIKELY`; fewer yields `NO_EVIDENCE_OF_MIXED_INFECTION`. This is
a screening heuristic, not a formal complexity-of-infection estimate; DEploid
(Zhu, Almagro-Garcia & McVean 2018) is available as a guarded wrapper for
formal deconvolution.

### 9.3 Relatedness

Pairwise identity-by-state (IBS) sharing is computed from shared alleles at
callable sites. Identical genotypes give `IBS2 = 1.0`; disjoint genotypes give
`IBS2 = 0.0`. A pairwise matrix is produced. hmmIBD (Schaffner et al. 2018) is
available as a guarded wrapper for formal identity-by-descent inference.

### 9.4 Phylogenetics

A neighbour-joining tree (Saitou & Nei 1987) is built from the pairwise distance
matrix and serialised in Newick format. IQ-TREE 2 (Minh et al. 2020) is
available as a guarded wrapper for maximum-likelihood inference with bootstrap
support (Hoang et al. 2018).

---

## 10. Quality control

QC metrics are defined in a machine-readable ontology
(`config/schemas/qc_ontology.json`, 17 metrics). Each metric declares its unit,
stage, direction (`min`, `max`, `range`, `none`), and whether its threshold is
provisional or validated. Evaluation returns `PASS`, `WARN`, `FAIL`, or `INFO`.
The run-level summary aggregates per-sample statuses and reports the overall
status as the weakest component. Every QC record carries the policy identifier,
so a QC verdict is traceable to the exact configuration that produced it.

---

## 11. Benchmarking and ablation

The benchmark harness simulates reads at controlled coverage and allele
fractions, calls them, and scores against the known truth. It reports SNP and
INDEL precision, recall, and F1, minority sensitivity across allele fraction and
depth, false positives in negative controls, calibration, and runtime/memory.
Confidence intervals are **Wilson score intervals** (Wilson 1927), which are
appropriate for proportions near 0 or 1 and avoid the degenerate intervals of
the Wald method (Brown, Cai & DasGupta 2001).

The benchmark harness can run simulations with known synthetic truth and
reports performance for those simulated conditions; this is a software/test
demonstration, not empirical biological validation. Some benchmark labels and
historical reports describe an ablation ladder, but those labels must not be read
as a description of the caller's inferential path: production calling uses the
exact one-sided count tail and complete-scan Bonferroni bound described in §5.3.
The repository does not bundle external truth data or a validated assay LoD.
Any comparisons to BCFtools, FreeBayes, or GATK are method-comparison context or
user-run benchmark options, not results established by the tests/demo described
here.

---

## 12. Reproducibility

Every run records a configuration hash (SHA-256 over the canonicalised
configuration), a parameter-registry hash, and a policy identifier of the form
`PVX-{mode}-{policy[:4]}-{config_hash[:12]}-{registry_hash[:8]}`. Writes are
atomic. Two runs with identical inputs and configuration produce byte-identical
evidence tables, which is verified by an integration test. Randomised
procedures (bootstrap, simulation) take an explicit seed.

---

## References

Brown, L.D., Cai, T.T. and DasGupta, A. (2001) 'Interval estimation for a
binomial proportion', *Statistical Science*, 16(2), pp. 101–133.

Chen, S., Zhou, Y., Chen, Y. and Gu, J. (2018) 'fastp: an ultra-fast all-in-one
FASTQ preprocessor', *Bioinformatics*, 34(17), pp. i884–i890.

Fisher, R.A. (1935) 'The logic of inductive inference', *Journal of the Royal
Statistical Society*, 98(1), pp. 39–82.

Fu, G.K., Hu, J., Wang, P.H. and Fodor, S.P.A. (2011) 'Counting individual DNA
molecules by the stochastic attachment of diverse labels', *Proceedings of the
National Academy of Sciences*, 108(22), pp. 9026–9031.

Gardner, M.J., Hall, N., Fung, E., White, O., Berriman, M., Hyman, R.W.,
Carlton, J.M., Pain, A., Nelson, K.E., Bowman, S. *et al.* (2002) 'Genome
sequence of the human malaria parasite *Plasmodium falciparum*', *Nature*,
419(6906), pp. 498–511.

Hoang, D.T., Chernomor, O., von Haeseler, A., Minh, B.Q. and Vinh, L.S. (2018)
'UFBoot2: improving the ultrafast bootstrap approximation', *Molecular Biology
and Evolution*, 35(2), pp. 518–522.

Jeffreys, H. (1946) 'An invariant form for the prior probability in estimation
problems', *Proceedings of the Royal Society A*, 186(1007), pp. 453–461.

Kinde, I., Wu, J., Papadopoulos, N., Kinzler, K.W. and Vogelstein, B. (2011)
'Detection and quantification of rare mutations with massively parallel
sequencing', *Proceedings of the National Academy of Sciences*, 108(23),
pp. 9530–9535.

Kivioja, T., Vähärautio, A., Karlsson, K., Bonke, M., Enge, M., Linnarsson, S.
and Taipale, J. (2011) 'Counting absolute numbers of molecules using unique
molecular identifiers', *Nature Methods*, 9(1), pp. 72–74.

Li, H. (2011) 'A statistical framework for SNP calling, mutation discovery,
association mapping and population genetical parameter estimation from
sequencing data', *Bioinformatics*, 27(21), pp. 2987–2993.

Li, H. (2013) 'Aligning sequence reads, clone sequences and assembly contigs
with BWA-MEM', *arXiv preprint*, arXiv:1303.3997.

MalariaGEN (2023) 'Pf7: an open dataset of *Plasmodium falciparum* genome
variation in 20,000 worldwide samples', *Wellcome Open Research*, 8, 22.

Minh, B.Q., Schmidt, H.A., Chernomor, O., Schrempf, D., Woodhams, M.D., von
Haeseler, A. and Lanfear, R. (2020) 'IQ-TREE 2: new models and efficient methods
for phylogenetic inference in the genomic era', *Molecular Biology and
Evolution*, 37(5), pp. 1530–1534.

Nei, M. and Li, W.H. (1979) 'Mathematical model for studying genetic variation
in terms of restriction endonucleases', *Proceedings of the National Academy of
Sciences*, 76(10), pp. 5269–5273.

Newman, A.M., Lovejoy, A.F., Klass, D.M., Kurtz, D.M., Chabon, J.J., Scherer, F.,
Stehr, H., Liu, C.L., Bratman, S.V., Say, C. *et al.* (2016) 'Integrated digital
error suppression for improved detection of circulating tumor DNA', *Nature
Biotechnology*, 34(5), pp. 547–555.

Nielsen, R., Paul, J.S., Albrechtsen, A. and Song, Y.S. (2011) 'Genotype and SNP
calling from next-generation sequencing data', *Nature Reviews Genetics*, 12(6),
pp. 443–451.

National Institute of Standards and Technology (2015) *Secure Hash Standard
(SHS)*. Federal Information Processing Standards Publication 180-4.
Gaithersburg, MD: NIST.

Saitou, N. and Nei, M. (1987) 'The neighbor-joining method: a new method for
reconstructing phylogenetic trees', *Molecular Biology and Evolution*, 4(4),
pp. 406–425.

Schaffner, S.F., Taylor, A.R., Wong, W., Munoz, D.E. and Neafsey, D.E. (2018)
'hmmIBD: software to infer pairwise identity by descent between haploid
genotypes', *Malaria Journal*, 17, 196.

Schmitt, M.W., Kennedy, S.R., Salk, J.J., Fox, E.J., Hiatt, J.B. and Loeb, L.A.
(2012) 'Detection of ultra-rare mutations by next-generation sequencing',
*Proceedings of the National Academy of Sciences*, 109(36), pp. 14508–14513.

Self, S.G. and Liang, K.Y. (1987) 'Asymptotic properties of maximum likelihood
estimators and likelihood ratio tests under nonstandard conditions', *Journal of
the American Statistical Association*, 82(398), pp. 605–610.

Tajima, F. (1989) 'Statistical method for testing the neutral mutation
hypothesis by DNA polymorphism', *Genetics*, 123(3), pp. 585–595.

Vasimuddin, M., Misra, S., Li, H. and Aluru, S. (2019) 'Efficient
architecture-aware acceleration of BWA-MEM for multicore systems', in
*Proceedings of the IEEE International Parallel and Distributed Processing
Symposium (IPDPS)*, pp. 314–324.

Watterson, G.A. (1975) 'On the number of segregating sites in genetical models
without recombination', *Theoretical Population Biology*, 7(2), pp. 256–276.

Wilson, E.B. (1927) 'Probable inference, the law of succession, and statistical
inference', *Journal of the American Statistical Association*, 22(158),
pp. 209–212.

Zhu, S.J., Almagro-Garcia, J. and McVean, G. (2018) 'Deconvoluting multiple
infections in *Plasmodium falciparum* from high throughput sequencing data',
*Bioinformatics*, 34(17), pp. i647–i655.
