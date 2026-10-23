"""Shared pytest fixtures for PfVarX."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pfvarx.config.schema import PfVarXConfig  # noqa: E402


def make_config(**overrides) -> PfVarXConfig:
    """Build a valid PROVISIONAL configuration for tests."""
    data = {
        "run": {"mode": "wgs", "intent": "research", "policy_status": "PROVISIONAL_DEVELOPMENT"},
        "reference": {
            "assembly_id": "Pf3D7",
            "fasta": "reference.fa",
            "sha256": "REQUIRED_SHA256",
        },
        "inputs": {"samplesheet": "samples.tsv"},
    }
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(data.get(key), dict):
            data[key].update(value)
        else:
            data[key] = value
    return PfVarXConfig.model_validate(data)


@pytest.fixture
def cfg() -> PfVarXConfig:
    return make_config()


@pytest.fixture
def synthetic(tmp_path):
    """A small synthetic dataset with three truth variants."""
    import pysam

    from pfvarx.synthetic import Variant, simulate_dataset, write_reference

    ref_path = str(tmp_path / "reference.fa")
    write_reference(ref_path, "Pf3D7_01_v3", 1200, seed=3)
    with pysam.FastaFile(ref_path) as fa:
        seq = fa.fetch("Pf3D7_01_v3").upper()

    def alt(b):
        return {"A": "T", "T": "A", "C": "G", "G": "C"}[b]

    variants = [
        Variant(300, seq[299], alt(seq[299]), 0.50),
        Variant(600, seq[599], alt(seq[599]), 0.20),
        Variant(900, seq[899], alt(seq[899]), 0.10),
    ]
    ds = simulate_dataset(
        str(tmp_path / "data"), length=1200, coverage=100, variants=variants, seed=3
    )
    return ds
