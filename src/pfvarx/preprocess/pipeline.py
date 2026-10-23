"""Command planning for preprocessing and alignment.

The planner produces an explicit, inspectable command plan for each sample
without executing anything.  This makes the workflow testable on machines that
do not have the external binaries installed, and makes every invocation
reproducible from the recorded plan.
"""

from __future__ import annotations

from pfvarx.config.schema import PfVarXConfig
from pfvarx.preprocess.tools import CommandPlan


def plan_fastq_route(
    cfg: PfVarXConfig,
    *,
    sample_id: str,
    r1: str,
    r2: str,
    out_dir: str,
    reference_index_prefix: str,
    threads: int | None = None,
) -> CommandPlan:
    """Plan fastp -> BWA-MEM2 -> samtools sort/index -> optional MarkDuplicates."""
    threads = threads or cfg.resources.threads_per_sample
    plan = CommandPlan(sample_id=sample_id)
    trimmed1 = f"{out_dir}/{sample_id}.trim.R1.fastq.gz"
    trimmed2 = f"{out_dir}/{sample_id}.trim.R2.fastq.gz"
    html = f"{out_dir}/{sample_id}.fastp.html"
    json = f"{out_dir}/{sample_id}.fastp.json"

    plan.add(
        "read_qc_trim",
        [
            "fastp", "-i", r1, "-I", r2, "-o", trimmed1, "-O", trimmed2,
            "--qualified_quality_phred", str(cfg.preprocessing.fastp_min_phred),
            "--length_required", str(cfg.preprocessing.fastp_min_length),
            "--thread", str(threads), "--html", html, "--json", json,
        ],
        [trimmed1, trimmed2, html, json],
        "fastp",
    )

    sorted_bam = f"{out_dir}/{sample_id}.sorted.bam"
    plan.add(
        "align",
        [
            "bwa-mem2", "mem", "-t", str(threads),
            "-R", f"@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA",
            reference_index_prefix, trimmed1, trimmed2,
        ],
        [f"{out_dir}/{sample_id}.sam"],
        "bwa-mem2",
    )
    plan.add(
        "sort_index",
        ["samtools", "sort", "-@", str(threads), "-o", sorted_bam,
         f"{out_dir}/{sample_id}.sam"],
        [sorted_bam, sorted_bam + ".bai"],
        "samtools",
    )

    if cfg.preprocessing.wgs_mark_duplicates and cfg.run.mode == "wgs":
        dedup = f"{out_dir}/{sample_id}.dedup.bam"
        metrics = f"{out_dir}/{sample_id}.duplicate_metrics.txt"
        plan.add(
            "mark_duplicates",
            ["picard", "MarkDuplicates", f"I={sorted_bam}", f"O={dedup}",
             f"M={metrics}", "CREATE_INDEX=true", "VALIDATION_STRINGENCY=LENIENT"],
            [dedup, dedup + ".bai", metrics],
            "picard",
        )
    return plan


def plan_umi_route(
    cfg: PfVarXConfig,
    *,
    sample_id: str,
    r1: str,
    r2: str,
    out_dir: str,
    reference_index_prefix: str,
    threads: int | None = None,
) -> CommandPlan:
    """Plan the UMI-aware branch: trim -> align -> UMI extract -> consensus -> realign.

    ``fgbio ExtractUmisFromBam`` operates on *aligned* reads, not raw FASTQ, so
    the reads are aligned first and the UMI tags are attached to the BAM.  The
    consensus reads are then re-aligned to produce the BAM consumed by the
    caller.
    """
    threads = threads or cfg.resources.threads_per_sample
    plan = CommandPlan(sample_id=sample_id)
    trimmed1 = f"{out_dir}/{sample_id}.trim.R1.fastq.gz"
    trimmed2 = f"{out_dir}/{sample_id}.trim.R2.fastq.gz"
    html = f"{out_dir}/{sample_id}.fastp.html"
    json = f"{out_dir}/{sample_id}.fastp.json"

    plan.add(
        "read_qc_trim",
        [
            "fastp", "-i", r1, "-I", r2, "-o", trimmed1, "-O", trimmed2,
            "--qualified_quality_phred", str(cfg.preprocessing.fastp_min_phred),
            "--length_required", str(cfg.preprocessing.fastp_min_length),
            "--thread", str(threads), "--html", html, "--json", json,
        ],
        [trimmed1, trimmed2, html, json],
        "fastp",
    )

    sorted_bam = f"{out_dir}/{sample_id}.sorted.bam"
    plan.add(
        "align",
        [
            "bwa-mem2", "mem", "-t", str(threads),
            "-R", f"@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA",
            reference_index_prefix, trimmed1, trimmed2,
        ],
        [f"{out_dir}/{sample_id}.sam"],
        "bwa-mem2",
    )
    plan.add(
        "sort_index",
        ["samtools", "sort", "-@", str(threads), "-o", sorted_bam,
         f"{out_dir}/{sample_id}.sam"],
        [sorted_bam, sorted_bam + ".bai"],
        "samtools",
    )

    umi_bam = f"{out_dir}/{sample_id}.umi.bam"
    plan.add(
        "umi_extract",
        ["fgbio", "ExtractUmisFromBam", f"--input={sorted_bam}",
         f"--output={umi_bam}",
         f"--read-structure={cfg.preprocessing.umi_design or 'M12S+T'}",
         "--single-tag=RX", "--molecular-index-tags=RX"],
        [umi_bam],
        "fgbio",
    )

    grouped = f"{out_dir}/{sample_id}.grouped.bam"
    plan.add(
        "group_consensus",
        ["fgbio", "GroupReadsByUmi", f"--input={umi_bam}",
         f"--output={grouped}", "--strategy=Adjacency", "--edits=1"],
        [grouped],
        "fgbio",
    )

    consensus_unmapped = f"{out_dir}/{sample_id}.consensus.unmapped.bam"
    plan.add(
        "consensus",
        ["fgbio", "CallMolecularConsensusReads", f"--input={grouped}",
         f"--output={consensus_unmapped}",
         f"--min-reads={cfg.preprocessing.umi_min_family_size}"],
        [consensus_unmapped],
        "fgbio",
    )

    consensus_sorted = f"{out_dir}/{sample_id}.consensus.sorted.bam"
    consensus_fastq = f"{out_dir}/{sample_id}.consensus.fastq.gz"
    plan.add(
        "consensus_to_fastq",
        ["samtools", "fastq", "-T", "RX", "-0", consensus_fastq,
         consensus_unmapped],
        [consensus_fastq],
        "samtools",
    )
    plan.add(
        "realign_consensus",
        [
            "bwa-mem2", "mem", "-t", str(threads),
            "-R", f"@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA",
            reference_index_prefix, consensus_fastq,
        ],
        [f"{out_dir}/{sample_id}.consensus.sam"],
        "bwa-mem2",
    )
    plan.add(
        "sort_consensus",
        ["samtools", "sort", "-@", str(threads), "-o", consensus_sorted,
         f"{out_dir}/{sample_id}.consensus.sam"],
        [consensus_sorted, consensus_sorted + ".bai"],
        "samtools",
    )
    return plan


def plan_route_for_sample(cfg: PfVarXConfig, route: str, **kwargs) -> CommandPlan:
    if route == "fastq_alignment":
        return plan_fastq_route(cfg, **kwargs)
    if route == "umi_consensus":
        return plan_umi_route(cfg, **kwargs)
    raise ValueError(f"No preprocessing plan for route {route!r}")
