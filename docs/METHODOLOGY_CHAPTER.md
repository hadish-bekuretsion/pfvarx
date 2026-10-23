# Chapter 3: Methodology

## 3.1 Introduction and methodological overview

This chapter sets out the methodological framework through which the present
study will characterise single-nucleotide polymorphisms (SNPs) and small
insertions and deletions (INDELs) in the *Plasmodium falciparum* genome, and
through which it will evaluate the reliability of that characterisation. The
research sits at the intersection of statistical genomics, computational
biology, and malaria epidemiology. Its central methodological claim is that the
detection of minority variants in polyclonal infections requires an explicit,
auditable statistical model rather than a threshold applied to an allele
fraction, and that the credibility of any such model depends on a validation
strategy that separates software correctness from biological accuracy.

The chapter is organised to mirror the analytical pipeline. It begins with the
research design and the epistemological commitments that shape it (§3.2), then
describes the data sources and sampling frame (§3.3), the reference and
annotation resources (§3.4), and the preprocessing and quality-control
procedures (§3.5). It then presents the variant-calling framework in detail
(§3.6), including the four count models and the decision-theoretic policy that
converts model output into calls (§3.7). Subsequent sections address variant
filtering and quality control (§3.8), minority-variant and polyclonal analysis
(§3.9), population-genomic and transmission analyses (§3.10–§3.12), functional
annotation (§3.13), computational reproducibility (§3.14), quality assurance
(§3.15), and the validation strategy (§3.16). The chapter closes with ethical
considerations (§3.17) and a discussion of limitations (§3.18).

A recurring theme is the distinction between *implementation* and *validation*.
A pipeline can be correctly implemented — its code can do what it claims — while
remaining scientifically unvalidated, because no empirical truth set has
demonstrated that its calls are correct on real data. This distinction, drawn
from software-engineering practice and applied here to genomic analysis, governs
the language of the entire chapter and the reporting conventions of the study.

## 3.2 Research design and epistemological stance

The software is intended for quantitative secondary analysis of existing
Illumina short-read data from *P. falciparum* isolates. This chapter distinguishes
the implemented pipeline from a proposed scientific study: this repository does
not document a completed analysis of characterized real-world samples. Simulations
and re-analysis of characterized datasets could investigate estimator behaviour,
but synthetic truth alone is not biological validation.

The epistemological stance is **critical realism** as applied to computational
measurement: a variant call is treated as an inference about an unobserved
biological state (the parasite genotype), mediated by a measurement process
(sequencing) and a model (the count distribution and decision rule). The
validity of the inference depends on the adequacy of the model to the
measurement process, not on the sophistication of the software. The implemented
caller uses configured fixed null parameters, not a fully calibrated distribution
that accounts for their estimation uncertainty or context-specific errors.

Three methodological commitments follow. First, **evidence before inference**:
every call must be traceable to a persisted record of the observations on which
it rests. Second, **statistical honesty**: a likelihood ratio must never be
presented as a posterior probability, and an uncalibrated statistic must never
be encoded as a Phred-scaled quality score. Third, **fail-closed behaviour**:
where an assumption cannot be met, the analysis must refuse to produce a call
rather than silently substitute a weaker method.

## 3.3 Data sources and sampling frame

The software is designed for *P. falciparum* whole-genome and targeted short-read
data, including datasets such as the MalariaGEN Pf7 resource (MalariaGEN, 2023).
That resource is contextual background, not a dataset analyzed or validated by
the tests/demo documented here. No external truth set, orthogonally characterized
dilution series, or validated assay LoD is bundled with this repository. Synthetic
data can exercise the implementation but cannot establish biological accuracy.

The sampling frame is defined by three inclusion criteria: (i) Illumina
paired-end short-read data; (ii) a demonstrable link to a pinned reference
assembly; and (iii) sufficient metadata to assign the sample to a sequencing
mode (whole-genome, targeted without unique molecular identifiers, or targeted
with unique molecular identifiers). Samples failing any criterion are excluded
and the exclusion is recorded. The unit of analysis is the sample for
within-host analyses and the variant site for population-genomic analyses.

These are proposed analysis criteria, not a record of samples analyzed in this
repository. No characterized external cohort or independent truth set is bundled.
Sample size, simulation replicates, and confidence intervals must be reported for
any future empirical or simulation study; the shipped tests/demo are not such a
study.

## 3.4 Reference genome and annotation resources

All analyses are anchored to a pinned *P. falciparum* reference assembly. The
reference genome was sequenced and annotated by the Malaria Genome Project
(Gardner *et al.*, 2002) and is distributed through PlasmoDB. The study will use
the Pf3D7 assembly (chromosomes `Pf3D7_01_v3` through `Pf3D7_14_v3`, together
with the apicoplast and mitochondrial genomes), consistent with the assembly
used by the MalariaGEN Pf7 resource (MalariaGEN, 2023).

The identity of the reference is recorded by both an assembly identifier and a
cryptographic checksum. The checksum convention follows the Secure Hash Standard
(NIST, 2015) and the reference-identifier practice of the GA4GH refget
specification; recording it prevents silent reference drift between runs, a
well-documented source of irreproducibility in genomic pipelines (Li, 2014).

Per-call annotation uses the internal engine. The configuration schema permits
`internal` or `auto` and rejects SnpEff selection; the SnpEff round-trip is not
implemented in the per-call path. With the necessary GFF3 and FASTA inputs, the
internal engine derives consequences; without them annotation may be unavailable.
Unit tests verify software behaviour, not scientific accuracy of consequence
predictions.

A **callable-region mask** is applied to exclude genomic regions where short-read
mapping is unreliable. These include the subtelomeric repeats and the
high-copy, highly polymorphic gene families (*var*, *rifin*, *stevor*) that are
known to be poorly resolved by short-read assembly and mapping (Gardner *et al.*,
2002; Otto *et al.*, 2018). Masking such regions is standard practice in
*P. falciparum* genomics and prevents spurious calls driven by multi-mapping
rather than biological variation.

## 3.5 Preprocessing and read-level quality control

### 3.5.1 Trimming and filtering

Adapter sequences and low-quality read tails are removed with fastp (Chen *et
al.*, 2018). The default thresholds are a minimum Phred score of 20, a minimum
retained length of 50 base pairs, and a maximum ambiguous-base fraction of 0.05.
These values follow conventional Illumina practice and are treated as
provisional, configurable parameters rather than fixed truths. The trimming step
is recorded so that the number of reads and bases removed is available for
quality control.

### 3.5.2 Alignment

Reads are aligned to the pinned reference with BWA-MEM2 (Li, 2013; Vasimuddin
*et al.*, 2019), which implements the BWA-MEM seed-and-extend algorithm with a
SIMD-accelerated backend. Mapping quality follows the BWA-MEM definition; reads
below a minimum mapping quality of 30 are excluded from evidence. The MAPQ 30
convention corresponds to the "uniquely mapped" threshold and is used by the
MalariaGEN Pf7 pipeline (MalariaGEN, 2023). Alignment statistics are summarised
with SAMtools (Li *et al.*, 2009; Danecek *et al.*, 2021).

### 3.5.3 Duplicate handling

For whole-genome data, polymerase-chain-reaction (PCR) duplicates are marked
with Picard MarkDuplicates and excluded from evidence. Duplicate exclusion is
important because PCR duplicates inflate apparent depth and can bias allele
fractions; the effect is well documented in Illumina library preparation (Aird
*et al.*, 2011). For targeted sequencing without unique molecular identifiers,
duplicate marking is disabled by default, because the fragment distribution is
not dominated by PCR in the same way; this is a provisional choice and is
reported as such.

### 3.5.4 Unique molecular identifiers

For targeted sequencing with unique molecular identifiers (UMIs), read families
are collapsed into consensus molecules with fgbio. UMI consensus reduces both
PCR and sequencing errors by requiring independent molecules to agree before a
base is accepted (Kivioja *et al.*, 2011; Fu *et al.*, 2011; Kinde *et al.*,
2011; Schmitt *et al.*, 2012; Newman *et al.*, 2016). A UMI design identifier is
mandatory for this mode; the pipeline fails closed if it is absent. The minimum
family size defaults to two molecules. UMI consensus is the principal mechanism
by which the study can detect variants below the error floor of conventional
sequencing, and its assumptions (that UMIs are unique and correctly attached)
are stated explicitly because they are not verifiable from the data alone.

## 3.6 The variant-calling framework

### 3.6.1 Scope and assumptions

The caller targets SNPs and small INDELs of 1–50 base pairs. It assumes a
haploid nuclear genome, consistent with the reference assembly (Gardner *et al.*,
2002), and it accommodates polyclonal infections by modelling the alternate
allele fraction as a continuous quantity rather than a binary genotype. It
contains no machine-learning component; every decision is a transparent
statistical rule. Structural variants, copy-number variation, large INDELs, and
non-Illumina platforms are outside scope.

### 3.6.2 Evidence extraction

Evidence extraction performs a CIGAR-aware pileup at each candidate site. It
records reference, alternate, and other base counts; informative depth; the
allele fraction with a Jeffreys credible interval (Jeffreys, 1946); mean base
quality; forward and reverse strand counts; mean read position and a
read-position bias statistic; a Fisher exact test for strand bias (Fisher,
1935); homopolymer run length; and left-alignment status for INDELs.

Three design choices warrant justification. First, overlapping mate pairs are
counted once, because a single DNA fragment observed twice is not two
independent observations. Second, low-mapping-quality reads, low-base-quality
bases, secondary alignments, and duplicate fragments are excluded **and counted
separately**, so that the evidence table records why depth was lost; this makes
the filtering auditable rather than opaque. Third, INDELs are represented in
anchored, left-aligned form, so that the same event is not counted twice under
different representations — a common source of false INDEL calls (Li, 2014).

The Jeffreys interval is used in preference to the Wald interval because the
latter behaves poorly near the boundaries `f = 0` and `f = 1`, which are exactly
the regions of interest for minority-variant detection (Brown, Cai and DasGupta,
2001). The Fisher exact test for strand bias follows standard contingency-table
logic and is appropriate for the small counts typical of low-fraction variants.

### 3.6.3 Count models

Four models are implemented, all evaluated in log space to avoid underflow at
high depth.

**Model A — binomial.** Under the null, alternate observations arise at
configured error rate `e`. Under H1, a true allele fraction `f` at or above the
configured boundary implies read-level ALT probability `f + (1-f)e`. The
binomial model assumes independent, identically distributed observations and
therefore understates variance when over-dispersion is present.

**Model B — beta-binomial.** The alternate-observation probability is drawn from
a Beta distribution with mean `mu` and over-dispersion `rho`, parameterised as
`alpha = mu(1 - rho)/rho` and `beta = (1 - mu)(1 - rho)/rho`. As `rho → 0` the
model converges to the binomial. The beta-binomial is an option for modelling
over-dispersion, not a universally validated default. Its `rho` requires
assay-specific calibration; the caller treats the configured value as fixed and
does not propagate estimation uncertainty. A literature value is not a
calibration for a particular assay.

**Model C — likelihood comparison and caller test.**
`count_models.evaluate_hypotheses` compares the error-only null with a composite
H1: true AF is at least `statistics.min_alt_fraction_h1`, and maximizes the
alternative likelihood over that allowed AF range. It reports
`2(logL1 - logL0)`. However, the caller does not use an LRT calibration as its
inferential p-value. `calling.policy` uses the exact one-sided upper tail of the
ALT count under the configured binomial or beta-binomial null; this value is
stored in the legacy field `lrt_pvalue`. Bootstrap and boundary-mixture LRT
helpers remain in the model module (the analytic boundary-mixture approximation
follows Self and Liang, 1987) but are not the caller's inference path.

For multiplicity, the caller applies Bonferroni using an upper bound over the
complete requested reference scan before candidate/evidence filtering. For each
anchor, the family bound includes three SNP alleles plus all possible anchored
insertion and deletion sequences through `candidate.max_small_indel_bp`. This
bound can be extremely conservative, especially because the possible small-
indel sequence space grows rapidly. The configuration fields
`exploratory_fdr_q` and `sensitivity_fdr_q` are legacy-named cutoffs applied to
the Bonferroni-adjusted p-value; this is not Benjamini–Hochberg/FDR control.
`bh_fdr` remains a helper, not the caller's procedure.
The correction is per sample and does not claim cohort-wide family-wise control
across multiple samples.

The count-tail p-value is exact conditional on the configured fixed null
parameters and informative depth. Error-rate and `rho` estimation uncertainty,
and context-specific error variation, are not represented by the caller.

**Model D — optional exploratory Bayesian quantity.** It is disabled by default.
When enabled, it uses the configured prior and a point AF boundary; it does not
integrate over a composite-AF prior, propagate null-parameter uncertainty, or
provide assay-validated posterior probability of biological truth.

The deliberate refusal to relabel a likelihood ratio as a posterior, or to
encode an uncalibrated statistic as a Phred quality score, is a methodological
safeguard. It prevents the common failure mode in which a test statistic is
presented as if it were a calibrated confidence, a practice that undermines
downstream interpretation.

## 3.7 Decision theory and the calling policy

Model output is converted into one of five decision categories by an explicit
policy: `NO_CALL`, `REFERENCE_SUPPORTED`, `ALT_DETECTED_EXPLORATORY`,
`ALT_PASS_VALIDATED`, and `ALT_UNCERTAIN`. The policy distinguishes **hard
gates**, which force `NO_CALL` (insufficient effective depth, a site outside the
callable region, or a reference-allele mismatch), from **review filters**, which
downgrade a call to `ALT_UNCERTAIN` (strand bias, homopolymer context, ambiguous
mapping, an invalid UMI model, or an uncalibrated model).

A `VALIDATED` configuration status, pinned reference, or error profile does not
by itself validate an assay. The current caller consumes no verifiable assay-
validation/LoD provenance and explicitly prevents `ALT_PASS_VALIDATED`; it cannot
produce a validated PASS. Calls can be exploratory, uncertain, reference-
supported, or no-call according to evidence and filters. This is a software
safeguard, not evidence that a biological error rate or decision threshold has
been validated.

The allele fraction is additionally classified into reporting bins: below 0.01
(below threshold), 0.01–0.05 (low-fraction exploratory), 0.05–0.20 (minority
candidate), 0.20–0.80 (intermediate fraction), and at or above 0.80
(predominant alternate). The intermediate bin is the principal signal for
polyclonal infection (§3.9).

## 3.8 Variant filtering and quality control

Two layers of quality control are applied. At the **site level**, the review
filters described in §3.7 flag calls for scrutiny, and the evidence table records
the reason for each flag. At the **sample and run level**, a machine-readable
quality-control ontology defines metrics with units, stages, directions (minimum,
maximum, range, or none), and thresholds. Each metric is evaluated as `PASS`,
`WARN`, `FAIL`, or `INFO`, and the run-level summary reports the weakest status
across components.

The separation of site-level and run-level QC is deliberate. Site-level filters
address the reliability of an individual call; run-level metrics address the
suitability of a sample or run for analysis. Conflating the two would allow a
globally poor run to pass because individual sites happened to meet thresholds,
or a good run to fail because a few sites were flagged.

## 3.9 Minority variants and polyclonal infection

*P. falciparum* infections are frequently polyclonal, so a single sample may
contain multiple genotypes and the alternate allele fraction at a site may lie
anywhere between zero and one (Manske *et al.*, 2012; Auburn *et al.*, 2011).
The study therefore treats the allele fraction as a continuous quantity and uses
the beta-binomial model (§3.6.3) to accommodate the resulting over-dispersion.

Samples are screened for polyclonality by counting variants in the intermediate
allele-fraction bin (0.20–0.80). Three or more such variants yields a
`MIXED_INFECTION_LIKELY` classification; fewer yields
`NO_EVIDENCE_OF_MIXED_INFECTION`. This is explicitly a screening heuristic, not
a formal complexity-of-infection estimate. Where formal deconvolution is
required, the study will apply DEploid (Zhu, Almagro-Garcia and McVean, 2018),
which estimates the number of constituent haplotypes and their proportions from
read counts. The heuristic and the formal method are reported separately so that
the reader can distinguish a screening flag from a modelled estimate.

## 3.10 Population-genomic analyses

Within-sample and between-sample diversity are quantified using standard
estimators. Nucleotide diversity `π` is the average pairwise difference per site
(Nei and Li, 1979). Watterson's theta `θ_W` is estimated from the number of
segregating sites (Watterson, 1975). Tajima's `D` compares the two and provides a
test of neutrality (Tajima, 1989). All estimators are computed over callable
bases, with missingness thresholds taken from the parameter registry, and are
reported with the number of callable bases so that the denominator is explicit.

Population structure and differentiation will be examined using allele-frequency
summaries and, where appropriate, model-based or distance-based clustering. The
study will report the number of variants and samples contributing to each
estimate, because population-genomic estimators are sensitive to missingness and
to the inclusion of low-quality sites.

## 3.11 Relatedness and transmission analysis

Pairwise relatedness is quantified by identity-by-state (IBS) sharing at callable
sites. Identical genotypes yield an IBS2 score of 1.0; disjoint genotypes yield
0.0. A pairwise matrix is produced for the cohort. Because IBS sharing does not
distinguish recent shared ancestry from chance allele sharing, the study will,
where the data permit, apply hmmIBD (Schaffner *et al.*, 2018), a hidden Markov
model that infers pairwise identity-by-descent between haploid genotypes. The
distinction between IBS and IBD is methodologically important: IBS is an
observation, whereas IBD is an inference about genealogy, and the two should not
be conflated in transmission analyses.

## 3.12 Phylogenetic analysis

A neighbour-joining tree (Saitou and Nei, 1987) is constructed from the pairwise
distance matrix and serialised in Newick format. Neighbour-joining is used as a
fast, deterministic baseline. Where maximum-likelihood inference is required,
the study will apply IQ-TREE 2 (Minh *et al.*, 2020) with ultrafast bootstrap
support (Hoang *et al.*, 2018). Trees are reported with their support values and
the substitution model used, and the study will state explicitly that a tree
built from a small number of variable sites has limited resolution.

## 3.13 Functional annotation and resistance markers

Per-call consequence annotation uses the internal translation engine; SnpEff
selection is rejected by the configuration schema and per-call SnpEff is not
implemented. Resistance markers are matched against the configured curated
catalogue. A match is an annotation from that catalogue, not a clinical
interpretation or independent validation of a resistance phenotype.

## 3.14 Computational reproducibility

Reproducibility is treated as a methodological requirement, not an
implementation detail. Every run records a configuration hash computed with
SHA-256 over the canonicalised configuration, a parameter-registry hash, and a
policy identifier that encodes both. Writes are atomic, so an interrupted run
leaves no truncated artefacts. Randomised procedures take explicit seeds. Two
runs with identical inputs and configuration produce byte-identical evidence
tables, a property verified by an automated test.

The analytical pipeline is orchestrated with Snakemake, which records the
directed acyclic graph of steps and supports resumption after interruption. The
Pixi manifest declares dependency version ranges, but the repository currently
has no committed lockfile; installs are not yet locked to a byte-identical
dependency solution. Container definitions are provided for Docker and
Apptainer/Singularity. These practices follow the reproducibility principles
articulated for computational biology by
Sandve *et al.* (2013) and the FAIR data principles (Wilkinson *et al.*, 2016).

## 3.15 Quality assurance procedures

Quality assurance operates at four levels. At the **code level**, the pipeline is
covered by unit and integration tests, and linting enforces a consistent style.
At the **data level**, inputs are validated for type, integrity, and completeness
before any compute is spent, and checksums are recorded. At the **analysis
level**, site-level filters and run-level QC metrics are applied as described in
§3.8. At the **reporting level**, every output carries a policy identifier and a
validation-status banner, so that a result cannot be presented without its
provenance and its validation status.

The four levels are complementary. Code tests cannot detect a biologically wrong
model; data validation cannot detect a coding error; analysis-level QC cannot
detect a misconfigured run; and reporting safeguards cannot detect a silently
degraded method. Only the combination provides defence in depth.

## 3.16 Validation strategy

Validation must distinguish software correctness from scientific accuracy.

**Software tests and demo.** Unit/integration tests check configured code paths,
interfaces, numerical routines, policy behaviour, and output contracts. The
synthetic demo/benchmark can score calls against its generated truth under the
specified simulation. Passing tests or scoring well on synthetic truth does not
demonstrate biological accuracy, performance on real samples, assay validation,
or a clinically relevant limit of detection.

**Scientific validation is not established in this repository.** Empirical
validation would require mode- and context-relevant characterized controls,
independent truth, pre-specified sensitivity/specificity criteria, and documented
assay-validation/LoD provenance. Calibration estimates of error rate and rho are
not sufficient on their own: this caller treats configured values as fixed and
does not propagate their estimation uncertainty or represent context-specific
errors. Comparisons to BCFtools, FreeBayes, or GATK are possible as method-
comparison work but are not results of the tests/demo described here. The
configured benchmark/ablation machinery does not change the production caller's
exact count-tail plus complete-scan Bonferroni inference path. Any proposed
external benchmarking, independent replication, and LoD study is future work,
not completed validation.

## 3.17 Ethical considerations

The study analyses secondary genomic data from *P. falciparum*, a parasite, and
therefore raises no human-subject consent issues in the usual sense. Nevertheless,
parasite genomic data are frequently linked to human clinical metadata, and the
study will handle such linkages in accordance with the data-governance policies
of the originating studies and the authors' institutions. Data will be stored and
processed within institutional infrastructure; the pipeline does not transmit data
externally. Where results could inform antimalarial treatment policy, the study
will present them as research findings and will not make clinical recommendations.

## 3.18 Limitations and mitigation

Important limitations are explicit. The caller has no verifiable assay-
validation/LoD provenance and cannot emit a validated PASS, even when the
configuration is labelled `VALIDATED`. The exact null-tail test conditions on
configured fixed error-rate and rho values; estimation uncertainty and context-
specific error rates are not represented. Bonferroni uses a complete pre-filter
scan upper bound that counts all possible configured small-indel alleles, which
can be extremely conservative and reduce power. Bayesian output is exploratory
and uses a point AF boundary. Synthetic tests/demo are not biological validation.
Other implemented summaries, including mixed-infection/relatedness, are screening
or descriptive outputs and should not be confused with formal estimates. The
internal annotation path is not a substitute for validated functional or clinical
interpretation.

## 3.19 Summary of methods

The software provides input routing, preprocessing, evidence extraction,
binomial/beta-binomial likelihoods, a decision policy, internal annotation,
downstream summaries, QC, and reporting. In the caller, `evaluate_hypotheses`
uses a composite H1 with AF at or above a configured boundary, while inferential
p-values are exact one-sided ALT-count tails under fixed configured null
parameters and corrected by Bonferroni over a conservative upper bound for the
complete pre-filter scan. The legacy `lrt_pvalue` / `fdr_q` names do not describe
LRT p-values or BH/FDR inference. Optional Bayesian output is exploratory and
point-boundary only. The caller cannot produce validated PASS without verifiable
assay-validation/LoD provenance. Tests and synthetic demo results demonstrate
software behaviour, not biological validation; empirical validation remains
future work.

---

## References

Aird, D., Ross, M.G., Chen, W.-S., Danielsson, M., Fennell, T., Russ, C., Jaffe,
D.B., Nusbaum, C. and Gnirke, A. (2011) 'Analyzing and minimizing PCR
amplification bias in Illumina sequencing libraries', *Genome Biology*, 12(2),
R18.

Auburn, S., Campino, S., Clark, T.G., Djimde, A.A., Zongo, I., Pinzon, M.,
Manske, M., MacInnis, B., Rockett, K.A., Wellems, T.E. and Kwiatkowski, D.P.
(2011) 'Characterization of within-host *Plasmodium falciparum* diversity using
next-generation sequence data', *PLoS ONE*, 6(9), e23691.

Brown, L.D., Cai, T.T. and DasGupta, A. (2001) 'Interval estimation for a
binomial proportion', *Statistical Science*, 16(2), pp. 101–133.

Chen, S., Zhou, Y., Chen, Y. and Gu, J. (2018) 'fastp: an ultra-fast all-in-one
FASTQ preprocessor', *Bioinformatics*, 34(17), pp. i884–i890.

Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O.,
Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M. and Li, H. (2021) 'Twelve
years of SAMtools and BCFtools', *GigaScience*, 10(2), giab008.

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

Li, H. (2014) 'Toward better understanding of artifacts in variant calling from
high-coverage samples', *Bioinformatics*, 30(20), pp. 2843–2851.

Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G.,
Abecasis, G., Durbin, R. and 1000 Genome Project Data Processing Subgroup (2009)
'The Sequence Alignment/Map format and SAMtools', *Bioinformatics*, 25(16),
pp. 2078–2079.

MalariaGEN (2023) 'Pf7: an open dataset of *Plasmodium falciparum* genome
variation in 20,000 worldwide samples', *Wellcome Open Research*, 8, 22.

Manske, M., Miotto, O., Campino, S., Auburn, S., Almagro-Garcia, J., Maslen, G.,
O'Brien, J., Djimde, A., Doumbo, O., Zongo, I. *et al.* (2012) 'Analysis of
*Plasmodium falciparum* diversity in natural infections by deep sequencing',
*Nature*, 487(7407), pp. 375–379.

Minh, B.Q., Schmidt, H.A., Chernomor, O., Schrempf, D., Woodhams, M.D., von
Haeseler, A. and Lanfear, R. (2020) 'IQ-TREE 2: new models and efficient methods
for phylogenetic inference in the genomic era', *Molecular Biology and
Evolution*, 37(5), pp. 1530–1534.

Miotto, O., Amato, R., Ashley, E.A., MacInnis, B., Almagro-Garcia, J., Amaratunga,
C., Lim, P., Mead, D., Oyola, S.O., Dhorda, M. *et al.* (2015) 'Genetic
architecture of artemisinin-resistant *Plasmodium falciparum*', *Nature Genetics*,
47(3), pp. 226–234.

Nei, M. and Li, W.H. (1979) 'Mathematical model for studying genetic variation in
terms of restriction endonucleases', *Proceedings of the National Academy of
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

Otto, T.D., Böhme, U., Sanders, M., Reid, A., Bruske, E.I., Duffy, C.W., Bull,
P.C., Pearson, R.D., Abdi, A., Dimonte, S. *et al.* (2018) 'Long read assemblies
of geographically dispersed *Plasmodium falciparum* isolates reveal highly
structured subtelomeres', *Wellcome Open Research*, 3, 52.

Saitou, N. and Nei, M. (1987) 'The neighbor-joining method: a new method for
reconstructing phylogenetic trees', *Molecular Biology and Evolution*, 4(4),
pp. 406–425.

Sandve, G.K., Nekrutenko, A., Taylor, J. and Hovig, E. (2013) 'Ten simple rules
for reproducible computational research', *PLoS Computational Biology*, 9(10),
e1003285.

Schaffner, S.F., Taylor, A.R., Wong, W., Munoz, D.E. and Neafsey, D.E. (2018)
'hmmIBD: software to infer pairwise identity by descent between haploid
genotypes', *Malaria Journal*, 17, 196.

Schmitt, M.W., Kennedy, S.R., Salk, J.J., Fox, E.J., Hiatt, J.B. and Loeb, L.A.
(2012) 'Detection of ultra-rare mutations by next-generation sequencing',
*Proceedings of the National Academy of Sciences*, 109(36), pp. 14508–14513.

Self, S.G. and Liang, K.Y. (1987) 'Asymptotic properties of maximum likelihood
estimators and likelihood ratio tests under nonstandard conditions', *Journal of
the American Statistical Association*, 82(398), pp. 605–610.

Tajima, F. (1989) 'Statistical method for testing the neutral mutation hypothesis
by DNA polymorphism', *Genetics*, 123(3), pp. 585–595.

Vasimuddin, M., Misra, S., Li, H. and Aluru, S. (2019) 'Efficient
architecture-aware acceleration of BWA-MEM for multicore systems', in
*Proceedings of the IEEE International Parallel and Distributed Processing
Symposium (IPDPS)*, pp. 314–324.

Watterson, G.A. (1975) 'On the number of segregating sites in genetical models
without recombination', *Theoretical Population Biology*, 7(2), pp. 256–276.

Wilkinson, M.D., Dumontier, M., Aalbersberg, I.J., Appleton, G., Axton, M.,
Baak, A., Blomberg, N., Boiten, J.-W., da Silva Santos, L.B., Bourne, P.E. *et
al.* (2016) 'The FAIR Guiding Principles for scientific data management and
stewardship', *Scientific Data*, 3, 160018.

Wilson, E.B. (1927) 'Probable inference, the law of succession, and statistical
inference', *Journal of the American Statistical Association*, 22(158),
pp. 209–212.

Zhu, S.J., Almagro-Garcia, J. and McVean, G. (2018) 'Deconvoluting multiple
infections in *Plasmodium falciparum* from high throughput sequencing data',
*Bioinformatics*, 34(17), pp. i647–i655.
