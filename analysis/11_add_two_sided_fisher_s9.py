#!/usr/bin/env python3
"""Add two-sided Fisher sensitivity columns to Repository Table S9.

The lineage-enrichment screen in the manuscript uses a one-sided Fisher exact
alternative because the analysis asks whether each lineage is over-represented
in the combined RTR Top100 context. This script keeps that primary screen
unchanged and adds two-sided Fisher exact p-values plus BH-FDR-adjusted values
as a sensitivity complement.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
S9_PATH = ROOT / "tables" / "S9_RTR_top100_lineage_enrichment_fisher_bh.csv"

TOTAL_N = 578
TOP_N = 100


def bh_fdr(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values, preserving input order."""
    p = [float(x) for x in pvalues]
    m = len(p)
    order = sorted(range(m), key=lambda i: p[i])
    adjusted_by_rank = [0.0] * m
    running_min = 1.0
    for rank_from_end, i in enumerate(reversed(order), start=1):
        rank = m - rank_from_end + 1
        adjusted = min(running_min, p[i] * m / rank)
        running_min = adjusted
        adjusted_by_rank[rank - 1] = adjusted
    out = [0.0] * m
    for rank, i in enumerate(order):
        out[i] = min(adjusted_by_rank[rank], 1.0)
    return out


def main() -> None:
    df = pd.read_csv(S9_PATH)
    required = {"lineage", "n_in_lineage", "n_in_top100", "fisher_p_one_sided"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"S9 is missing required columns: {sorted(missing)}")

    fisher_p_two_sided: list[float] = []
    fisher_p_one_sided_recomputed: list[float] = []

    for _, row in df.iterrows():
        lineage = row["lineage"]
        a = int(row["n_in_top100"])
        n_lineage = int(row["n_in_lineage"])
        if not (0 <= a <= min(TOP_N, n_lineage)):
            raise ValueError(f"Invalid Top100/lineage counts for {lineage!r}")

        c = n_lineage - a
        b = TOP_N - a
        d = (TOTAL_N - TOP_N) - c
        if min(b, c, d) < 0:
            raise ValueError(f"Invalid 2x2 table for {lineage!r}: {[[a, c], [b, d]]}")

        # Rows are lineage / not-lineage; columns are in Top100 / not-in-Top100.
        table = [[a, c], [b, d]]
        _, p_one = fisher_exact(table, alternative="greater")
        _, p_two = fisher_exact(table, alternative="two-sided")
        fisher_p_one_sided_recomputed.append(float(p_one))
        fisher_p_two_sided.append(float(p_two))

    max_abs_diff = max(
        abs(float(old) - float(new))
        for old, new in zip(df["fisher_p_one_sided"], fisher_p_one_sided_recomputed)
    )
    if max_abs_diff > 1e-12:
        raise RuntimeError(
            "Existing one-sided Fisher values were not reproduced; "
            f"maximum absolute difference = {max_abs_diff:.3g}"
        )

    df["fisher_p_two_sided"] = fisher_p_two_sided
    df["BH_FDR_two_sided"] = bh_fdr(fisher_p_two_sided)
    df.to_csv(S9_PATH, index=False)
    print(f"[written] {S9_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
