"""RTR concordance and Primary-only enrichment analyses.

Recomputes two combined-RTR robustness checks from outputs/tables/cell_line_RTR_scores_with_quadrants_and_variants.csv: Secondary-vs-Primary rank concordance across dual-evaluable cell lines, and Primary-only over-representation in the combined RTR Top100.
"""
import math
import numpy as np
import pandas as pd

DF = pd.read_csv("outputs/tables/cell_line_RTR_scores_with_quadrants_and_variants.csv")
sec = DF["secondary_RTR_delta_any_repair_minus_removal"]
pri = DF["primary_RTR_delta_any_repair_minus_removal"]
dual = DF[sec.notna() & pri.notna()]
prionly = DF[sec.isna() & pri.notna()]

# Spearman correlation (Pearson on percentile ranks)
rx = dual["secondary_RTR_delta_any_repair_minus_removal"].rank().values
ry = dual["primary_RTR_delta_any_repair_minus_removal"].rank().values
rho = float(np.corrcoef(rx, ry)[0, 1])
print(f"[RTR concordance] Spearman rho = {rho:.3f} (n={len(dual)})")

# Primary-only enrichment: hypergeometric / Fisher exact (exact, math.comb)
N, K, n = len(DF), len(prionly), 100
top = DF.sort_values("RTR_score_combined_rank_0_1", ascending=False).head(n)
k = int((top["secondary_RTR_delta_any_repair_minus_removal"].isna()
         & top["primary_RTR_delta_any_repair_minus_removal"].notna()).sum())
C = math.comb
pmf = lambda x: C(K, x) * C(N - K, n - x) / C(N, n)
xmax = min(K, n)
p_ge = sum(pmf(x) for x in range(k, xmax + 1))
p0 = pmf(k)
p_two = sum(pmf(x) for x in range(xmax + 1) if pmf(x) <= p0 * (1 + 1e-9))
a, b, c, d = k, n - k, K - k, (N - n) - (K - k)
OR = (a * d) / (b * c)
print(f"[Primary-only enrichment] Top100 primary-only {k}/{n}; background {K}/{N}; "
      f"P(X>={k})={p_ge:.3g}; Fisher2={p_two:.3g}; OR={OR:.2f}")


# Persist outputs for repository-table synchronization when run as a script.
from pathlib import Path as _Path
outdir = _Path("outputs/tables"); outdir.mkdir(parents=True, exist_ok=True)
pd.DataFrame([dict(
    metric="spearman_rho_secondary_vs_primary_RTR_contrast",
    value=round(rho, 3),
    n=len(dual),
    note="Spearman correlation between the Secondary and Primary per-dataset repair-minus-removal percentile ranks across dual-evaluable cell lines; computed from outputs/tables/cell_line_RTR_scores_with_quadrants_and_variants.csv. Spearman is invariant to monotone transform, so this equals the correlation of the per-dataset percentile ranks that are averaged in the combined RTR score. Used as a rank-averaging concordance check for the combined RTR score.",
)]).to_csv(outdir/"RTR_secondary_primary_rank_concordance.csv", index=False)
pd.DataFrame([
    dict(quantity="top100_primary_only_count", value=str(k), note="Primary-only cell lines among the combined RTR Top100"),
    dict(quantity="top100_size", value=str(n), note="combined RTR Top100 size"),
    dict(quantity="top100_primary_only_fraction", value=f"{k/n:.2f}", note=f"{k}/{n}"),
    dict(quantity="background_primary_only_count", value=str(K), note="Primary-only cell lines among all evaluable combined-score cell lines"),
    dict(quantity="background_total", value=str(N), note="evaluable combined-score cell lines"),
    dict(quantity="background_primary_only_fraction", value=f"{K/N:.4f}", note=f"{K}/{N} = {100*K/N:.1f}%"),
    dict(quantity="hypergeometric_one_sided_P_X_ge_31", value=f"{p_ge:.1e}", note="P(X>=31) under hypergeometric null (over-representation)"),
    dict(quantity="fisher_exact_two_sided_p", value=f"{p_two:.1e}", note="two-sided Fisher exact on 2x2 [[31,69],[67,411]]"),
    dict(quantity="odds_ratio", value=f"{OR:.2f}", note="odds ratio for Primary-only enrichment in Top100"),
    dict(quantity="note_overall", value="Computed from outputs/tables/cell_line_RTR_scores_with_quadrants_and_variants.csv. Primary-only = secondary contrast NaN and primary contrast present. Used to quantify Primary-only over-representation in the combined Top100. Mitigated by the dual-evaluable-only sensitivity check (69/100 overlap; Bone+Soft Tissue aggregate remained enriched).", note=""),
]).to_csv(outdir/"RTR_top100_primary_only_enrichment.csv", index=False)
# keep S16 deliverable copies synchronized when the script is run from the repository root
s16 = _Path("tables/S16_RTR_combined_score_robustness"); s16.mkdir(parents=True, exist_ok=True)
for _name in ["RTR_secondary_primary_rank_concordance.csv", "RTR_top100_primary_only_enrichment.csv"]:
    (s16 / _name).write_text((outdir / _name).read_text(encoding="utf-8"), encoding="utf-8")
