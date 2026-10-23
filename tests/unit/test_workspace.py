"""Tests for three-phase workspace planning."""

from __future__ import annotations

import pytest

from pfvarx.errors import ConfigError
from pfvarx.workspace import (
    DIR_ALIGNMENTS,
    DIR_OBJECTIVES,
    create,
    create_input_directories,
    plan_input_directories,
    plan_layout,
    prepare,
    validate,
)


def test_fastq_route_creates_alignment_directories(tmp_path):
    layout = plan_layout(tmp_path, objectives=["diversity"], routes=["fastq_alignment"])
    assert DIR_ALIGNMENTS in layout.directories
    assert f"{DIR_OBJECTIVES}/diversity" in layout.directories


def test_bam_route_does_not_create_alignment_directories(tmp_path):
    layout = plan_layout(tmp_path, objectives=["diversity"], routes=["alignment_reuse"])
    assert DIR_ALIGNMENTS not in layout.directories
    assert f"{DIR_OBJECTIVES}/diversity" in layout.directories


def test_vcf_route_does_not_create_alignment_directories(tmp_path):
    layout = plan_layout(tmp_path, objectives=["coi"], routes=["variant_import"])
    assert DIR_ALIGNMENTS not in layout.directories


def test_unselected_objectives_get_no_directory(tmp_path):
    layout = plan_layout(tmp_path, objectives=["diversity"], routes=["alignment_reuse"])
    objective_dirs = {
        d.rsplit("/", 1)[-1]
        for d in layout.directories
        if d.startswith(DIR_OBJECTIVES + "/")
    }
    assert objective_dirs == {"diversity"}
    for absent in ("relatedness", "phylogenetics", "drug_resistance", "population_structure"):
        assert absent not in objective_dirs


def test_every_directory_has_a_rationale(tmp_path):
    layout = plan_layout(tmp_path, objectives=["diversity"], routes=["alignment_reuse"])
    for directory in layout.directories:
        assert layout.rationale.get(directory), f"{directory} has no rationale"


def test_planning_is_pure_and_creates_nothing(tmp_path):
    plan_layout(tmp_path, objectives=["diversity"], routes=["alignment_reuse"])
    assert list(tmp_path.iterdir()) == []


def test_validate_warns_about_existing_directories_without_side_effects(tmp_path):
    logs = tmp_path / "logs"
    logs.mkdir()
    sentinel = logs / "keep.txt"
    sentinel.write_text("existing output")
    layout = plan_layout(tmp_path, objectives=["diversity"], routes=["alignment_reuse"])
    warnings = validate(layout, overwrite=False)
    assert warnings and "overwrite is false" in warnings[0]
    assert sentinel.read_text() == "existing output"
    assert sorted(path.name for path in tmp_path.iterdir()) == ["logs"]
    assert not (tmp_path / "config").exists()


def test_overwrite_validation_does_not_clear_existing_output(tmp_path):
    logs = tmp_path / "logs"
    logs.mkdir()
    sentinel = logs / "keep.txt"
    sentinel.write_text("existing output")
    layout = plan_layout(tmp_path, objectives=["diversity"], routes=["alignment_reuse"])
    assert validate(layout, overwrite=True) == []
    create(layout)
    assert sentinel.read_text() == "existing output"


def test_validate_is_silent_for_a_fresh_root(tmp_path):
    layout = plan_layout(tmp_path, objectives=["diversity"], routes=["alignment_reuse"])
    assert validate(layout, overwrite=False) == []


def test_validate_rejects_a_file_at_the_output_root(tmp_path):
    target = tmp_path / "results"
    target.write_text("not a directory")
    layout = plan_layout(target, objectives=["diversity"], routes=["alignment_reuse"])
    with pytest.raises(ConfigError):
        validate(layout, overwrite=False)


def test_create_makes_only_planned_directories(tmp_path):
    layout = plan_layout(tmp_path, objectives=["diversity"], routes=["alignment_reuse"])
    created = create(layout)
    for directory in layout.directories:
        assert (tmp_path / directory).is_dir()
    assert not (tmp_path / DIR_ALIGNMENTS).exists()
    assert len(created) == len(layout.directories)


def test_prepare_returns_created_directories(tmp_path):
    layout, created, warnings = prepare(
        tmp_path, objectives=["coi"], routes=["variant_import"], overwrite=False
    )
    assert created
    assert warnings == []
    assert (tmp_path / DIR_OBJECTIVES / "coi").is_dir()


def test_prepare_is_idempotent_for_existing_directories(tmp_path):
    prepare(tmp_path, objectives=["coi"], routes=["variant_import"], overwrite=False)
    _, created, warnings = prepare(
        tmp_path, objectives=["coi"], routes=["variant_import"], overwrite=False
    )
    assert created == []
    assert warnings


def test_as_dict_round_trips(tmp_path):
    layout = plan_layout(tmp_path, objectives=["diversity"], routes=["alignment_reuse"])
    payload = layout.as_dict()
    assert payload["root"] == str(tmp_path)
    assert payload["directories"] == list(layout.directories)
    assert isinstance(payload["rationale"], dict)


def test_input_directory_setup_creates_no_outputs(tmp_path):
    from tests.conftest import make_config

    data = tmp_path / "custom-inputs"
    results = tmp_path / "custom-results"
    logs = tmp_path / "custom-logs"
    work = tmp_path / "custom-work"
    cfg = make_config(
        reference={"fasta": str(data / "reference" / "Pf3D7.fa")},
        inputs={
            "root": str(data),
            "raw_dir": str(data / "raw"),
            "prepared_dir": str(data / "prepared"),
            "metadata_dir": str(data / "metadata"),
            "samplesheet": str(data / "prepared" / "samples.tsv"),
        },
        outputs={
            "root": str(results),
            "logs_dir": str(logs),
            "work_dir": str(work),
        },
    )

    planned = plan_input_directories(cfg)
    created = create_input_directories(cfg)

    assert set(created) == set(planned)
    assert (data / "raw").is_dir()
    assert (data / "prepared").is_dir()
    assert (data / "metadata").is_dir()
    assert (data / "reference").is_dir()
    assert not results.exists()
    assert not logs.exists()
    assert not work.exists()
    assert create_input_directories(cfg) == []


def test_input_directory_setup_rejects_output_overlap(tmp_path):
    from tests.conftest import make_config

    output = tmp_path / "results"
    cfg = make_config(
        inputs={"root": str(tmp_path), "raw_dir": str(output)},
        outputs={"root": str(output)},
    )
    with pytest.raises(ConfigError, match="overlaps a configured output"):
        create_input_directories(cfg)
    assert not output.exists()
