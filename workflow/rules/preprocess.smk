# Stage 1 — preprocessing and alignment (FASTQ route only).
#
# fastp -> BWA-MEM2 -> samtools sort -> optional Picard MarkDuplicates.
# For TBS-with-UMI mode the UMI-aware branch (fgbio) is used instead.
# Existing BAM/CRAM inputs never reach these rules.

rule fastp_trim:
    input:
        r1=lambda wc: sample_record(wc.sample).fastq_1,
        r2=lambda wc: sample_record(wc.sample).fastq_2,
    output:
        r1=f"{WORK}/{{sample}}/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{WORK}/{{sample}}/trimmed/{{sample}}_R2.fastq.gz",
        json=f"{WORK}/{{sample}}/trimmed/{{sample}}.fastp.json",
        html=f"{WORK}/{{sample}}/trimmed/{{sample}}.fastp.html",
    params:
        min_phred=config.get("preprocessing", {}).get("fastp_min_phred", 20),
        min_len=config.get("preprocessing", {}).get("fastp_min_length", 50),
        adapter="--detect_adapter_for_pe"
        if config.get("preprocessing", {}).get("adapter_autodetect", True)
        else "",
    threads: THREADS
    resources:
        mem_mb=2048,
    log:
        f"{LOGS}/{{sample}}/fastp.log",
    shell:
        """
        mkdir -p {WORK}/{wildcards.sample}/trimmed {LOGS}/{wildcards.sample}
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              --qualified_quality_phred {params.min_phred} \
              --length_required {params.min_len} \
              {params.adapter} \
              --json {output.json} --html {output.html} --thread {threads} 2> {log}
        """


rule bwa_mem2_align:
    input:
        r1=f"{WORK}/{{sample}}/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{WORK}/{{sample}}/trimmed/{{sample}}_R2.fastq.gz",
    output:
        bam=f"{WORK}/{{sample}}/aligned/{{sample}}.sorted.bam",
        bai=f"{WORK}/{{sample}}/aligned/{{sample}}.sorted.bam.bai",
    params:
        ref=REF_FASTA,
        rg=lambda wc: f"-R '@RG\\tID:{wc.sample}\\tSM:{wc.sample}\\tPL:ILLUMINA'",
    threads: THREADS
    resources:
        mem_mb=MEM_MB,
    log:
        f"{LOGS}/{{sample}}/bwa.log",
    shell:
        """
        mkdir -p {WORK}/{wildcards.sample}/aligned
        bwa-mem2 mem -t {threads} {params.rg} {params.ref} {input.r1} {input.r2} \
          | samtools sort -@ {threads} -m {resources.mem_mb}M -o {output.bam} - 2> {log}
        samtools index {output.bam}
        """


rule mark_duplicates:
    input:
        bam=f"{WORK}/{{sample}}/aligned/{{sample}}.sorted.bam",
    output:
        bam=f"{WORK}/{{sample}}/aligned/{{sample}}.dedup.bam",
        bai=f"{WORK}/{{sample}}/aligned/{{sample}}.dedup.bam.bai",
        metrics=f"{WORK}/{{sample}}/aligned/{{sample}}.dup_metrics.txt",
    params:
        enabled=config.get("preprocessing", {}).get("wgs_mark_duplicates", True),
    threads: THREADS
    resources:
        mem_mb=MEM_MB,
    log:
        f"{LOGS}/{{sample}}/picard.log",
    shell:
        """
        if [ "{params.enabled}" = "True" ]; then
            picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} \
                REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT 2> {log}
        else
            cp {input.bam} {output.bam}
            printf 'mark_duplicates disabled by configuration\n' > {output.metrics}
        fi
        samtools index {output.bam}
        """


# ---------------------------------------------------------------------------
# UMI-aware branch (tbs_umi).  Only used when MODE == "tbs_umi" and the sample
# arrives as FASTQ.  fgbio operates on aligned reads, so the trimmed reads are
# aligned first; consensus reads are converted back to FASTQ and re-aligned.
# ---------------------------------------------------------------------------

rule umi_extract:
    input:
        bam=f"{WORK}/{{sample}}/aligned/{{sample}}.sorted.bam",
    output:
        bam=f"{WORK}/{{sample}}/aligned/{{sample}}.umi.bam",
    params:
        design=config.get("preprocessing", {}).get("umi_design", "M12S+T"),
    log:
        f"{LOGS}/{{sample}}/umi_extract.log",
    shell:
        """
        fgbio ExtractUmisFromBam --input={input.bam} --output={output.bam} \
            --read-structure={params.design} --single-tag=RX \
            --molecular-index-tags=RX 2> {log}
        """


rule group_consensus:
    input:
        bam=f"{WORK}/{{sample}}/aligned/{{sample}}.umi.bam",
    output:
        bam=f"{WORK}/{{sample}}/aligned/{{sample}}.grouped.bam",
    log:
        f"{LOGS}/{{sample}}/group_consensus.log",
    shell:
        """
        fgbio GroupReadsByUmi --input={input.bam} --output={output.bam} \
            --strategy=Adjacency --edits=1 2> {log}
        """


rule call_consensus:
    input:
        bam=f"{WORK}/{{sample}}/aligned/{{sample}}.grouped.bam",
    output:
        bam=f"{WORK}/{{sample}}/aligned/{{sample}}.consensus.unmapped.bam",
    params:
        min_reads=config.get("preprocessing", {}).get("umi_min_family_size", 2),
    log:
        f"{LOGS}/{{sample}}/call_consensus.log",
    shell:
        """
        fgbio CallMolecularConsensusReads --input={input.bam} \
            --output={output.bam} --min-reads={params.min_reads} 2> {log}
        """


rule consensus_to_fastq:
    input:
        bam=f"{WORK}/{{sample}}/aligned/{{sample}}.consensus.unmapped.bam",
    output:
        fastq=f"{WORK}/{{sample}}/aligned/{{sample}}.consensus.fastq.gz",
    log:
        f"{LOGS}/{{sample}}/consensus_to_fastq.log",
    shell:
        """
        samtools fastq -T RX -0 {output.fastq} {input.bam} 2> {log}
        """


rule realign_consensus:
    input:
        fastq=f"{WORK}/{{sample}}/aligned/{{sample}}.consensus.fastq.gz",
    output:
        sam=f"{WORK}/{{sample}}/aligned/{{sample}}.consensus.sam",
    params:
        ref=REF_FASTA,
        rg=lambda wc: f"-R '@RG\\tID:{wc.sample}\\tSM:{wc.sample}\\tPL:ILLUMINA'",
    threads: THREADS
    log:
        f"{LOGS}/{{sample}}/realign_consensus.log",
    shell:
        """
        bwa-mem2 mem -t {threads} {params.rg} {params.ref} {input.fastq} \
            > {output.sam} 2> {log}
        """


rule sort_consensus:
    input:
        sam=f"{WORK}/{{sample}}/aligned/{{sample}}.consensus.sam",
    output:
        bam=f"{WORK}/{{sample}}/aligned/{{sample}}.consensus.sorted.bam",
        bai=f"{WORK}/{{sample}}/aligned/{{sample}}.consensus.sorted.bam.bai",
    threads: THREADS
    log:
        f"{LOGS}/{{sample}}/sort_consensus.log",
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.sam} 2> {log}
        samtools index {output.bam}
        """
