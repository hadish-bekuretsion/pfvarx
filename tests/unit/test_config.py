"""Configuration, schema and parameter-registry tests."""

from __future__ import annotations

import pytest

from pfvarx.config.parameters import (
    get_parameter,
    registry_hash,
    registry_records,
    validate_parameter,
)
from pfvarx.config.schema import PfVarXConfig
from pfvarx.errors import ConfigError
from tests.conftest import make_config


def test_config_hash_is_stable_and_order_independent():
    a = make_config()
    b = make_config()
    assert a.config_hash() == b.config_hash()
    assert len(a.config_hash()) == 64


def test_policy_id_encodes_mode_and_hash():
    cfg = make_config()
    pid = cfg.policy_id()
    assert pid.startswith("PVX-wgs-")
    assert cfg.config_hash()[:12] in pid


def test_unknown_field_rejected():
    with pytest.raises((ConfigError, ValueError)):
        PfVarXConfig.model_validate(
            {
                "run": {"mode": "wgs"},
                "reference": {"assembly_id": "x", "fasta": "f", "sha256": "REQUIRED_SHA256"},
                "inputs": {"samplesheet": "s.tsv"},
                "bogus_section": {"a": 1},
            }
        )


def test_invalid_mode_rejected():
    with pytest.raises((ConfigError, ValueError)):
        make_config(run={"mode": "rna"})


def test_unimplemented_snpeff_engine_is_rejected_in_config():
    with pytest.raises((ConfigError, ValueError)):
        make_config(annotation={"engine": "snpeff"})


def test_tbs_umi_requires_umi_design():
    with pytest.raises((ConfigError, ValueError)):
        PfVarXConfig.model_validate(
            {
                "run": {"mode": "tbs_umi"},
                "reference": {"assembly_id": "x", "fasta": "f", "sha256": "REQUIRED_SHA256"},
                "inputs": {"samplesheet": "s.tsv"},
                "preprocessing": {"umi_design": None},
            }
        )


def test_validated_requires_pinned_assembly():
    with pytest.raises((ConfigError, ValueError)):
        PfVarXConfig.model_validate(
            {
                "run": {"mode": "wgs", "policy_status": "VALIDATED"},
                "reference": {"assembly_id": "REQUIRED_PINNED", "fasta": "f", "sha256": "REQUIRED_SHA256"},
                "inputs": {"samplesheet": "s.tsv"},
            }
        )


def test_bayesian_requires_posterior_threshold():
    with pytest.raises((ConfigError, ValueError)):
        make_config(statistics={"bayesian_enabled": True})


def test_reporting_fractions_must_be_ordered():
    with pytest.raises((ConfigError, ValueError)):
        make_config(reporting={"minority_candidate_fraction": 0.9, "intermediate_fraction": 0.2})


def test_parameter_registry_has_required_metadata():
    for rec in registry_records():
        for field in ("name", "unit", "datatype", "default", "rationale", "source", "policy_status"):
            assert field in rec, f"{rec.get('name')} missing {field}"


def test_parameter_validation_rejects_out_of_range():
    p = get_parameter("preprocessing.fastp_min_phred")
    with pytest.raises((ConfigError, ValueError)):
        validate_parameter("preprocessing.fastp_min_phred", p.maximum + 1)


def test_parameter_validation_accepts_default():
    p = get_parameter("preprocessing.fastp_min_phred")
    assert validate_parameter("preprocessing.fastp_min_phred", p.default) == p.default


def test_registry_hash_stable():
    assert registry_hash() == registry_hash()
    assert len(registry_hash()) == 64
