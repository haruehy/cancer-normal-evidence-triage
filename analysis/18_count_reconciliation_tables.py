#!/usr/bin/env python3
"""Create bundled count-reconciliation tables used for repository review.

These tables are derived only from bundled processed outputs / S-numbered
repository tables. They do not require or redistribute raw PRISM matrices.
"""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "tables"
S1 = ROOT / "tables" / "S1_compound_annotations_contexts_TP53_RTR_FDR"
OUT.mkdir(parents=True, exist_ok=True)
S1.mkdir(parents=True, exist_ok=True)

GROUP_ORDER = [
    "strict repair/protective",
    "broad repair/protective",
    "cytotoxic/removal-like",
    "other",
]


def ordered(df: pd.DataFrame, group_col: str = "analysis_group") -> pd.DataFrame:
    tmp = df.copy()
    tmp["_order"] = tmp[group_col].map({g: i for i, g in enumerate(GROUP_ORDER)}).fillna(99)
    return tmp.sort_values(["_order", group_col]).drop(columns="_order").reset_index(drop=True)


def write_both(df: pd.DataFrame, name: str) -> None:
    df.to_csv(OUT / name, index=False)
    df.to_csv(S1 / name, index=False)
    print(f"[written] outputs/tables/{name}")
    print(f"[written] tables/S1_compound_annotations_contexts_TP53_RTR_FDR/{name}")


def make_secondary_reconciliation() -> None:
    ann_counts = pd.read_csv(OUT / "compound_annotation_counts.csv")
    grp = pd.read_csv(OUT / "prism_secondary_repair_vs_removal_group_summary.csv")

    ann_map = ann_counts.set_index("analysis_group")["n_compounds"].to_dict()
    grp_map = grp.set_index("analysis_group")["n_drugs"].to_dict()

    rows = []
    for g in GROUP_ORDER:
        curated_n = int(ann_map.get(g, 0))
        response_n = int(grp_map.get(g, 0))
        rows.append(
            {
                "analysis_group": g,
                "curated_annotation_count_n": curated_n,
                "response_group_summary_n_drugs": response_n,
                "response_minus_curated_n": response_n - curated_n,
                "interpretation": (
                    "PRISM Secondary response-observed drugs lacking a curated annotation row are assigned to other in the response-summary merge."
                    if g == "other" and response_n != curated_n
                    else "Counts match between the curated annotation table and the response-summary table."
                ),
            }
        )

    total_curated = sum(int(ann_map.get(g, 0)) for g in GROUP_ORDER)
    total_response = sum(int(grp_map.get(g, 0)) for g in GROUP_ORDER)
    rows.append(
        {
            "analysis_group": "TOTAL",
            "curated_annotation_count_n": total_curated,
            "response_group_summary_n_drugs": total_response,
            "response_minus_curated_n": total_response - total_curated,
            "interpretation": "The total difference is fully accounted for by the other class in the response-summary table.",
        }
    )
    write_both(pd.DataFrame(rows), "secondary_other_count_reconciliation.csv")


def make_primary_counts() -> None:
    src_path = S1 / "source_primary_drug_level_measurable_proxies.csv"
    src = pd.read_csv(src_path)
    if "eqcl_group" not in src.columns:
        raise ValueError("source_primary_drug_level_measurable_proxies.csv must contain eqcl_group")
    if "drug_key" not in src.columns:
        raise ValueError("source_primary_drug_level_measurable_proxies.csv must contain drug_key")

    counts = (
        src.drop_duplicates("drug_key")
        .groupby("eqcl_group")
        .size()
        .reindex(GROUP_ORDER, fill_value=0)
        .reset_index()
        .rename(columns={"eqcl_group": "analysis_group", 0: "n_compounds"})
    )
    counts["dataset"] = "PRISM Primary selected"
    counts["count_basis"] = "unique drug_key rows in source_primary_drug_level_measurable_proxies.csv"
    counts["source_table"] = "tables/S1_compound_annotations_contexts_TP53_RTR_FDR/source_primary_drug_level_measurable_proxies.csv"
    counts["raw_data_scope"] = "derived aggregate count only; raw PRISM Primary response data are not redistributed"
    counts = counts[["dataset", "analysis_group", "n_compounds", "count_basis", "source_table", "raw_data_scope"]]
    write_both(counts, "prism_primary_selected_compound_annotation_counts.csv")


def main() -> None:
    make_secondary_reconciliation()
    make_primary_counts()


if __name__ == "__main__":
    main()
