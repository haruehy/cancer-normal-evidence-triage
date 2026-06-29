# S6 Table
AUC clipping, strict-only limiting-case summaries (see strict rows in AUC_clip_sensitivity_main_results), and RTR repair-preservation/removal-vulnerability weight sensitivity. broad_subclass_membership.csv lists the broad-group subclass assignments used for subclass-exclusion sensitivity; the per-exclusion median-AUC recomputation requires the raw PRISM secondary matrix (see data/README). Rows with blank broad_subclass values are Strict R/P compounds that also belong to the broad superset; the manuscript subclass counts refer to the 75 broad repair/protective rows.

Additional clipped-TP53 sensitivity support:
- `AUC_clip_TP53_cell_line_delta.csv` reports the clipped Secondary cell-line repair-minus-cytotoxic delta comparison for the `TP53 LoF/high-impact` versus `TP53 mutation-not-called` contrast. This file supports the manuscript's clipped-TP53 p-value statement (MW-AUC 0.5185, p=0.5389).
