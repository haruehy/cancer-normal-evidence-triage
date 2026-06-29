#!/usr/bin/env python3
"""Rule-based consistency audit for single-rater compound-class annotations.

This script does not assign manuscript classes and does not estimate inter-rater
reliability. It checks whether the released definition-compound annotations for
the Strict R/P, Broad R/P, and cytotoxic/removal-like classes contain at least
one transparent support term in the compound name, MOA, target, or annotation text.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ANN = ROOT / "data" / "curated_compound_annotations.csv"
OUT = ROOT / "outputs" / "tables"
AUDIT = ROOT / "tables" / "curation_consistency_audit"

STRICT_PATTERNS = [
    r"antioxidant", r"radical", r"osmolyte", r"cytoprotect", r"ergothioneine",
    r"taurine", r"trolox", r"benfotiamine", r"glutathione", r"mesna", r"chlorogenic",
]
BROAD_PATTERNS = STRICT_PATTERNS + [
    r"anti.?inflamm", r"cyclooxygenase", r"\bcox\b", r"lipoxygenase", r"\blox\b",
    r"nf.?kb", r"\bikk\b", r"tnf", r"tumor necrosis factor", r"glucocorticoid",
    r"corticosteroid", r"ppar", r"nitric oxide", r"nos", r"prostaglandin",
    r"histamine", r"jak", r"interleukin", r"leukotriene", r"cannabinoid", r"trpv", r"faah",
]
CYTOTOXIC_PATTERNS = [
    r"topoisomerase", r"\bparp\b", r"\batr\b", r"\bchk", r"\bwee", r"aurora",
    r"\bplk", r"\bkif", r"kinesin", r"tubulin", r"microtubule", r"proteasome",
    r"hsp90", r"bcl", r"mcl", r"iap", r"apoptosis", r"mitotic", r"cdk",
    r"hdac", r"mek", r"egfr", r"alkylat", r"anthracycline", r"tax", r"vinca",
    r"vinblastine", r"camptothecin", r"dna", r"kinase inhibitor", r"\bp53\b",
    r"\btp53\b", r"\bmdm\b", r"caspase", r"casp[0-9]",
]
PATTERNS = {
    "strict repair/protective": STRICT_PATTERNS,
    "broad repair/protective": BROAD_PATTERNS,
    "cytotoxic/removal-like": CYTOTOXIC_PATTERNS,
}


def matched_terms(text: str, patterns: list[str]) -> list[str]:
    text = str(text).lower()
    return [p for p in patterns if re.search(p, text)]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)
    ann = pd.read_csv(ANN)
    if "annotation_text" not in ann.columns:
        ann["annotation_text"] = ann[["name", "moa", "target"]].fillna("").astype(str).agg(" | ".join, axis=1)

    detail_rows = []
    summary_rows = []
    for group, patterns in PATTERNS.items():
        sub = ann[(ann["analysis_group"] == group) & (ann["is_definition_compound"].astype(bool))].copy()
        for _, row in sub.iterrows():
            search_text = " | ".join(str(row.get(c, "")) for c in ["name", "moa", "target", "annotation_text"])
            terms = matched_terms(search_text, patterns)
            detail_rows.append({
                "analysis_group": group,
                "drug_key": row.get("drug_key", ""),
                "name": row.get("name", ""),
                "moa": row.get("moa", ""),
                "target": row.get("target", ""),
                "support_term_found": bool(terms),
                "matched_support_terms": ";".join(terms),
            })
        passed = sum(1 for r in detail_rows if r["analysis_group"] == group and r["support_term_found"])
        total = len(sub)
        summary_rows.append({
            "analysis_group": group,
            "definition_compounds_checked": total,
            "definition_compounds_with_support_term": passed,
            "definition_compounds_without_support_term": total - passed,
            "support_term_coverage": passed / total if total else 0.0,
            "audit_interpretation": "support-term consistency check only; not inter-rater reliability",
        })

    detail = pd.DataFrame(detail_rows)
    summary = pd.DataFrame(summary_rows)
    for path in [OUT / "curation_consistency_audit_definition_compounds.csv", AUDIT / "curation_consistency_audit_definition_compounds.csv"]:
        detail.to_csv(path, index=False)
    for path in [OUT / "curation_consistency_audit_summary.csv", AUDIT / "curation_consistency_audit_summary.csv"]:
        summary.to_csv(path, index=False)

    readme = AUDIT / "README.md"
    readme.write_text(
        "# Curation consistency audit\n\n"
        "This audit checks whether definition compounds in the released single-rater "
        "compound annotation table contain at least one transparent support term in "
        "the compound name, mechanism-of-action, target, or annotation text. It is "
        "a consistency/transparency check only. It is not an inter-rater reliability "
        "estimate, does not replace independent curation, and does not assign the "
        "manuscript analysis groups.\n\n"
        "Outputs:\n"
        "- `curation_consistency_audit_summary.csv`\n"
        "- `curation_consistency_audit_definition_compounds.csv`\n",
        encoding="utf-8",
    )
    print("Wrote curation consistency audit outputs")


if __name__ == "__main__":
    main()
