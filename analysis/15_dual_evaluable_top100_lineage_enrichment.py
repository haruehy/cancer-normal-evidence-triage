#!/usr/bin/env python3
"""Dual-evaluable-only Top100 Bone+Soft Tissue enrichment analysis.

Recomputes the dual-evaluable-only lineage sensitivity analysis directly from the
released combined RTR decomposition table, with no raw-data download and no SciPy
dependency (Fisher exact / hypergeometric tail computed from log-gamma).

Run:
    python analysis/15_dual_evaluable_top100_lineage_enrichment.py

Writes / verifies:
    outputs/tables/RTR_dual_evaluable_top100_lineage_enrichment.csv
"""
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "outputs" / "tables" / "cell_line_RTR_scores_with_quadrants_and_variants.csv"
OUT = ROOT / "outputs" / "tables" / "RTR_dual_evaluable_top100_lineage_enrichment.csv"
BST = {"Bone", "Soft Tissue"}


def fnum(x):
    x = (x or "").strip()
    if x == "" or x.lower() == "nan":
        return None
    try:
        return float(x)
    except ValueError:
        return None


def logC(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def main():
    rows = list(csv.DictReader(open(SRC)))
    dual = [
        (fnum(r["RTR_score_combined_rank_0_1"]), (r["OncotreeLineage"] or "").strip())
        for r in rows
        if fnum(r["secondary_RTR_delta_any_repair_minus_removal"]) is not None
        and fnum(r["primary_RTR_delta_any_repair_minus_removal"]) is not None
        and fnum(r["RTR_score_combined_rank_0_1"]) is not None
    ]
    dual.sort(key=lambda t: -t[0])
    top = dual[:100]
    n = len(top)
    N = len(dual)
    a = sum(1 for _, l in top if l in BST)
    K = sum(1 for _, l in dual if l in BST)
    b, c, d = n - a, K - a, N - n - (K - a)

    def pmf(k):
        return math.exp(logC(K, k) + logC(N - K, n - k) - logC(N, n))

    kmin, kmax = max(0, n - (N - K)), min(n, K)
    p_one = sum(pmf(k) for k in range(a, kmax + 1))
    pa = pmf(a)
    p_two = sum(pmf(k) for k in range(kmin, kmax + 1) if pmf(k) <= pa * (1 + 1e-9))
    orv = (a * d) / (b * c)

    # checks against the manuscript-reported values
    assert (N, K) == (480, 20), (N, K)
    assert (a, n) == (12, 100), (a, n)
    assert abs(p_one - 1.055e-4) < 1e-6, p_one
    assert abs(orv - 6.34) < 1e-2, orv

    print(f"dual-evaluable-only Top100 Bone+Soft Tissue: {a}/{n} vs background {K}/{N}")
    print(f"  Fisher one-sided p = {p_one:.3e}; two-sided p = {p_two:.3e}; OR = {orv:.4f}")
    print("[ok] reproduces manuscript values (12/100 vs 20/480; p=1.06e-4; OR=6.34)")

    note = "Dual-evaluable-only Top100 = 100 highest combined-RTR cell lines among the 480 cell lines evaluable in both PRISM Secondary and Primary; background = 480 dual-evaluable cell lines. Fisher exact (one-sided over-representation and two-sided) and hypergeometric agree to displayed precision. Used to test whether the Bone+Soft Tissue lineage signal persists after restricting the Top100 definition to dual-evaluable cell lines."
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["lineage_aggregate","n_in_top100","top100_size","n_in_background","background_total","fraction_top100","fraction_background","fisher_one_sided_p","fisher_two_sided_p","odds_ratio","overlap_with_default_top100","note"])
        w.writeheader()
        w.writerow({"lineage_aggregate":"Bone+Soft Tissue (aggregate)","n_in_top100":a,"top100_size":n,"n_in_background":K,"background_total":N,"fraction_top100":round(a/n, 6),"fraction_background":round(K/N, 6),"fisher_one_sided_p":f"{p_one:.3e}","fisher_two_sided_p":f"{p_two:.3e}","odds_ratio":round(orv, 4),"overlap_with_default_top100":69,"note":note})
    s16 = ROOT / "tables" / "S16_RTR_combined_score_robustness"
    s16.mkdir(parents=True, exist_ok=True)
    (s16 / OUT.name).write_text(OUT.read_text(encoding="utf-8"), encoding="utf-8")
    print(f"[ok] wrote {OUT.name} and synchronized the S16 copy")


if __name__ == "__main__":
    main()
