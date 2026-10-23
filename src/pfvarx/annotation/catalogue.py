"""Curated drug-resistance catalogue loader and matcher.

The catalogue is data, not code: it is versioned, carries citations and an
explicit ``clinical_use`` flag, and is loaded from YAML so that it can be
updated without changing the software.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from pathlib import Path

import yaml

from pfvarx.errors import ConfigError

DEFAULT_CATALOGUE = Path(__file__).with_name("pf_resistance.yaml")
_AA_SUBSTITUTION = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY*])([1-9][0-9]*)([ACDEFGHIKLMNPQRSTVWY*])$")


@dataclass
class ResistanceMarker:
    gene_id: str
    gene_name: str
    drug: str
    aa_change: str
    evidence_level: str
    notes: str
    references: list[str]

    def as_dict(self) -> dict:
        return {
            "gene_id": self.gene_id,
            "gene_name": self.gene_name,
            "drug": self.drug,
            "aa_change": self.aa_change,
            "evidence_level": self.evidence_level,
            "notes": self.notes,
            "references": "; ".join(self.references),
        }


@dataclass
class ResistanceCatalogue:
    version: str
    assembly: str
    clinical_use: str
    markers: list[ResistanceMarker]

    def by_gene(self, gene_id: str) -> list[ResistanceMarker]:
        return [m for m in self.markers if m.gene_id == gene_id]

    def match(self, gene_id: str, aa_change: str | None) -> list[ResistanceMarker]:
        if not aa_change:
            return []
        normalized = aa_change.removeprefix("p.")
        substitution = _AA_SUBSTITUTION.fullmatch(normalized)
        if not substitution:
            return []
        hits = []
        for m in self.by_gene(gene_id):
            marker_change = m.aa_change.removeprefix("p.")
            if marker_change == normalized:
                hits.append(m)
            elif (
                substitution.group(1) != substitution.group(3)
                and "nonsynonymous" in m.aa_change.casefold()
            ):
                # These catalogue descriptions lack explicit residue and (for
                # PfKelch13) domain constraints. Preserve them only as broad,
                # low-confidence review candidates; never inherit "established".
                hits.append(replace(m, evidence_level="low_confidence_candidate"))
        return hits


def load_catalogue(path: str | None = None) -> ResistanceCatalogue:
    p = Path(path) if path else DEFAULT_CATALOGUE
    try:
        with open(p) as fh:
            data = yaml.safe_load(fh)
    except FileNotFoundError as exc:
        raise ConfigError(
            f"Resistance catalogue not found: {p}",
            cause=str(exc),
            action="Provide a valid catalogue path or reinstall the package data.",
        ) from exc
    markers = [
        ResistanceMarker(
            gene_id=m["gene_id"],
            gene_name=m["gene_name"],
            drug=m["drug"],
            aa_change=str(m["aa_change"]),
            evidence_level=m.get("evidence_level", "candidate"),
            notes=m.get("notes", ""),
            references=m.get("references", []),
        )
        for m in data.get("markers", [])
    ]
    return ResistanceCatalogue(
        version=str(data.get("version", "unknown")),
        assembly=str(data.get("assembly", "unknown")),
        clinical_use=str(data.get("clinical_use", "NOT_FOR_CLINICAL_USE")),
        markers=markers,
    )
