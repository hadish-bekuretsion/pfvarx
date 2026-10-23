"""Variant annotation: internal consequence engine, SnpEff wrapper, catalogue."""

from pfvarx.annotation.annotator import AnnotatedSite, annotate_calls
from pfvarx.annotation.catalogue import ResistanceCatalogue, ResistanceMarker, load_catalogue
from pfvarx.annotation.consequence import Consequence, annotate_variant
from pfvarx.annotation.gff import CDS, GeneModel, parse_gff3, revcomp, translate

__all__ = [
    "AnnotatedSite", "CDS", "Consequence", "GeneModel", "ResistanceCatalogue",
    "ResistanceMarker", "annotate_calls", "annotate_variant", "load_catalogue",
    "parse_gff3", "revcomp", "translate",
]
