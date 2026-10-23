# Shared helpers for the PfVarX Snakemake workflow.

import os
from pathlib import Path

PFX = config.get("project", "PfVarX")
RESULTS = config.get("outputs", {}).get("root", "results")
LOGS = config.get("outputs", {}).get("logs_dir") or str(Path(RESULTS) / "logs")
WORK = config.get("outputs", {}).get("work_dir") or str(Path(RESULTS) / "intermediate")
THREADS = config.get("resources", {}).get("threads_per_sample", 2)
MEM_MB = config.get("resources", {}).get("mem_mb_per_heavy_rule", 6144)
REF_FASTA = config.get("reference", {}).get("fasta")
CALLABLE_BED = config.get("reference", {}).get("callable_mask_bed")


def sample_record(sample):
    return next(s for s in _samples if s.sample_id == sample)


def route_of(sample):
    return ROUTES[sample].route


def pfvarx_cli():
    """Return the pfvarx entry point, preferring the installed console script."""
    return config.get("pfvarx_cli", "pfvarx")


def q(x):
    """Shell-quote a path or value."""
    import shlex

    return shlex.quote(str(x))
