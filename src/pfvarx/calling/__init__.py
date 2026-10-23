"""Variant calling: evidence extraction, statistical policy and output."""

from pfvarx.calling.caller import CallingResult, call_sample, load_callable_mask
from pfvarx.calling.evidence import (
    Candidate,
    SiteEvidence,
    discover_indels,
    discover_snps,
    extract_site_evidence,
)
from pfvarx.calling.policy import (
    ALL_BINS,
    ALL_DECISIONS,
    ALL_FILTERS,
    DECISION_ALT_EXPLORATORY,
    DECISION_ALT_PASS,
    DECISION_ALT_UNCERTAIN,
    DECISION_NO_CALL,
    DECISION_REFERENCE_SUPPORTED,
    POLICY_VERSION,
    SiteCall,
    apply_fdr,
    evaluate_site,
    minority_bin,
    summarise_calls,
)
from pfvarx.calling.writer import EVIDENCE_COLUMNS, write_evidence_table, write_vcf

__all__ = [
    "ALL_BINS", "ALL_DECISIONS", "ALL_FILTERS",
    "DECISION_ALT_EXPLORATORY", "DECISION_ALT_PASS", "DECISION_ALT_UNCERTAIN",
    "DECISION_NO_CALL", "DECISION_REFERENCE_SUPPORTED", "POLICY_VERSION",
    "Candidate", "CallingResult", "EVIDENCE_COLUMNS", "SiteCall", "SiteEvidence",
    "apply_fdr", "call_sample", "discover_indels", "discover_snps",
    "evaluate_site", "extract_site_evidence", "load_callable_mask",
    "minority_bin", "summarise_calls", "write_evidence_table", "write_vcf",
]
