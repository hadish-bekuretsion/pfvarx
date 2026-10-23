"""Tests for the objective registry and dependency resolution."""

from __future__ import annotations

import pytest

from pfvarx.errors import ConfigError
from pfvarx.objectives import (
    ALL_OBJECTIVES,
    OBJECTIVES,
    format_table,
    normalise,
    resolve,
)


def test_registry_covers_every_documented_objective():
    assert set(OBJECTIVES) == set(ALL_OBJECTIVES)
    assert len(ALL_OBJECTIVES) == 8


def test_every_objective_declares_a_scope_and_stages():
    for obj_id, obj in OBJECTIVES.items():
        assert obj.id == obj_id
        assert obj.scope in {"per_sample", "run"}
        assert obj.requires, f"{obj_id} declares no required stages"
        assert obj.label and obj.summary


def test_aliases_resolve_to_canonical_ids():
    assert normalise("resistance") == "drug_resistance"
    assert normalise("mixed") == "mixed_infection"
    assert normalise("popstruct") == "population_structure"
    assert normalise("phylo") == "phylogenetics"
    assert normalise("spatial") == "spatial_temporal"
    assert normalise("Drug Resistance") == "drug_resistance"


def test_unknown_objective_raises_config_error_with_guidance():
    with pytest.raises(ConfigError) as exc:
        normalise("not_a_real_objective")
    message = str(exc.value)
    assert "Unknown analysis objective" in message
    assert "pfvarx objectives" in message


def test_resolution_is_order_independent_and_deterministic():
    a = resolve(["diversity", "coi"], have_catalogue=False)
    b = resolve(["coi", "diversity"], have_catalogue=False)
    assert a.selected == b.selected
    # Documented registry order, not user order.
    assert a.selected == [o for o in ALL_OBJECTIVES if o in {"diversity", "coi"}]


def test_resolution_unions_required_stages():
    resolved = resolve(["diversity", "mixed_infection"], have_catalogue=False)
    # Both objectives need the calling stage; the union is what will be executed.
    assert resolved.stages == ["call"]
    assert resolved.multi_sample is False
    assert resolved.needs_metadata is False


def test_drug_resistance_requires_catalogue():
    with pytest.raises(ConfigError) as exc:
        resolve(["drug_resistance"], have_catalogue=False)
    assert "catalogue" in str(exc.value).lower()

    resolved = resolve(["drug_resistance"], have_catalogue=True)
    assert resolved.needs_catalogue is True
    assert "annotate" in resolved.stages


def test_metadata_objective_requires_metadata():
    with pytest.raises(ConfigError) as exc:
        resolve(["population_structure"], have_metadata=False)
    assert "metadata" in str(exc.value).lower()

    resolved = resolve(["population_structure"], have_metadata=True)
    assert resolved.needs_metadata is True


def test_multi_sample_objective_requires_two_samples():
    with pytest.raises(ConfigError) as exc:
        resolve(["relatedness"], n_samples=1)
    assert "at least two samples" in str(exc.value).lower()

    resolved = resolve(["relatedness"], n_samples=3)
    assert resolved.multi_sample is True


def test_empty_selection_is_allowed_and_empty():
    resolved = resolve([], have_catalogue=True)
    assert resolved.is_empty
    assert resolved.selected == []
    assert resolved.stages == []


def test_duplicate_selection_is_collapsed():
    resolved = resolve(["diversity", "diversity", "coi"], have_catalogue=False)
    assert resolved.selected.count("diversity") == 1


def test_format_table_lists_every_objective_and_the_selector_hint():
    text = format_table()
    for obj_id in ALL_OBJECTIVES:
        assert obj_id in text
    assert "--objectives" in text
    assert "objectives:" in text


def test_as_dict_is_json_ready():
    resolved = resolve(["diversity"], have_catalogue=False)
    payload = resolved.as_dict()
    assert payload["selected"] == ["diversity"]
    assert isinstance(payload["labels"], list)
    assert isinstance(payload["stages"], list)
