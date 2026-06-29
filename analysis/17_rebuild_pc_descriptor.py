#!/usr/bin/env python3
"""Rebuild the final four-anchor P_C descriptor from reproducible PRISM-derived ranks.

The pair-priority source table contains both a binary repair-risk rank and a refined repair-risk rank. The final manuscript descriptor uses only the refined PRISM-derived survival-preservation rank. The binary rank is retained only as source-table provenance and is not used here.
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
S2_DIR = ROOT / "tables" / "S2_cancer_normal_evidence_table"
OUT = ROOT / "outputs" / "tables"

ANCHOR_DISPLAY = {
    "l-ergothioneine": "L-ergothioneine",
    "taurine": "Taurine",
    "trolox": "Trolox",
    "benfotiamine": "Benfotiamine",
}


def main() -> None:
    source_path = S2_DIR / "source_repair_candidates_used.csv"
    if not source_path.exists():
        raise FileNotFoundError(source_path)
    src = pd.read_csv(source_path)
    required = {"repair_drug_key", "repair_candidate", "refined_repair_risk_rank"}
    missing = required - set(src.columns)
    if missing:
        raise ValueError("source_repair_candidates_used.csv is missing required columns: " + ", ".join(sorted(missing)))
    src["_candidate_norm"] = src["repair_candidate"].astype(str).str.strip().str.lower()
    anchors = src[src["_candidate_norm"].isin(ANCHOR_DISPLAY)].copy()
    if len(anchors) != 4:
        raise ValueError(f"Expected four anchor R/P candidates, found {len(anchors)}")
    anchors["candidate"] = anchors["_candidate_norm"].map(ANCHOR_DISPLAY)
    anchors["final_PC_value"] = anchors["refined_repair_risk_rank"].astype(float)
    anchors["PC_table_value"] = anchors["final_PC_value"].round(3)
    anchors = anchors.sort_values("final_PC_value", ascending=False).reset_index(drop=True)
    anchors["PC_rank_within_four"] = [f"{i+1}/4" for i in range(len(anchors))]
    out = anchors[["candidate", "repair_drug_key", "refined_repair_risk_rank", "final_PC_value", "PC_table_value", "PC_rank_within_four"]].copy()
    out.insert(3, "PC_formula", "refined_repair_risk_rank only; binary_repair_risk_rank is retained only as provenance and is not used in final P_C")
    out["PC_generation_record"] = "rebuilt by analysis/17_rebuild_pc_descriptor.py from tables/S2_cancer_normal_evidence_table/source_repair_candidates_used.csv"
    OUT.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT / "final_pc_descriptor_rebuilt.csv", index=False)
    out.to_csv(S2_DIR / "pc_descriptor_rebuilt.csv", index=False)
    print("Rebuilt final P_C descriptor for four anchor R/P candidates")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
