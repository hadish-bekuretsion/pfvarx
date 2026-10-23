# Stage 2 — evidence extraction, statistical calling and variant QC.
#
# A single rule produces the evidence table and VCF for each sample, choosing
# the correct path from the sample's route:
#   * alignment_reuse -> call the supplied BAM/CRAM directly (no realignment)
#   * fastq_alignment -> call the BAM produced by the preprocessing rules
#   * variant_import  -> import a pre-called VCF/BCF/gVCF (no realignment)
#
# Dispatching inside one rule avoids two rules claiming the same output, which
# would otherwise let the import path shadow the caller for BAM inputs.

def _bam_for(sample):
    route = route_of(sample)
    if route == "alignment_reuse":
        return sample_record(sample).bam
    if MODE == "tbs_umi":
        return f"{WORK}/{sample}/aligned/{sample}.consensus.sorted.bam"
    return f"{WORK}/{sample}/aligned/{sample}.dedup.bam"


def _evidence_inputs(wc):
    """Return only the input the sample's route actually needs.

    A non-applicable input is returned as an empty list so Snakemake does not
    try to resolve it.  The `run:` block reads only the key relevant to the
    route.
    """
    if route_of(wc.sample) == "variant_import":
        return {"vcf": sample_record(wc.sample).vcf, "bam": []}
    return {"vcf": [], "bam": _bam_for(wc.sample)}


rule produce_evidence:
    input:
        plan=f"{RESULTS}/routing_plan.json",
        bam=lambda wc: _evidence_inputs(wc)["bam"],
        vcf=lambda wc: _evidence_inputs(wc)["vcf"],
    output:
        vcf=f"{RESULTS}/{{sample}}/{{sample}}.pfvarx.vcf",
        evidence=f"{RESULTS}/{{sample}}/{{sample}}.evidence.tsv",
    params:
        cli=pfvarx_cli(),
        cfg=CONFIG_FILE,
        ref=REF_FASTA,
        bed_arg=(f"--callable-bed {CALLABLE_BED}" if CALLABLE_BED else ""),
        sample=lambda wc: wc.sample,
    threads: THREADS
    resources:
        mem_mb=MEM_MB,
    log:
        f"{LOGS}/{{sample}}/call.log",
    run:
        import os

        os.makedirs(f"{LOGS}/{wildcards.sample}", exist_ok=True)
        if route_of(wildcards.sample) == "variant_import":
            shell(
                r"""
                {params.cli} import \
                    --config {params.cfg} \
                    --vcf {input.vcf} \
                    --sample-id {wildcards.sample} \
                    --out {RESULTS} 2> {log}
                """
            )
        else:
            shell(
                r"""
                {params.cli} call \
                    --config {params.cfg} \
                    --bam {input.bam} \
                    --sample-id {params.sample} \
                    --fasta {params.ref} \
                    {params.bed_arg} \
                    --out {RESULTS} 2> {log}
                """
            )
