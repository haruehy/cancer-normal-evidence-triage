# Changelog

## PRISM Secondary AUC reconciliation (post-v30)

### Fixed: stale compound classification in PRISM-Secondary-derived outputs
- `outputs/tables/prism_secondary_cell_line_context.csv` and related processed
  outputs had been generated under an older compound classification (strict
  repair/protective = 3, broad repair/protective = 142, cytotoxic/removal-like
  = 193) that predated the current curated annotation (strict = 6, broad = 75,
  cytotoxic = 157, matching Table 1). These outputs were never regenerated
  after the annotation was updated.
- `analysis/03_auc_clip_sensitivity.py` was also missing the
  `passed_str_profiling` filter applied in `analysis/02_prism_repair_vs_
  removal_analysis.py`; corrected to match.
- `01_build_compound_annotations.py`, `02_prism_repair_vs_removal_
  analysis.py`, and `03_auc_clip_sensitivity.py` were re-run in sequence
  against the current annotation and the raw `secondary-screen-dose-
  response-curve-parameters.csv`. All PRISM-Secondary-derived
  `outputs/tables/*` files now reflect the current (6/75/157) classification.

Manuscript values updated accordingly (Sections 4.3 and 4.7):

| Location | Old | Corrected |
|---|---|---|
| 4.3 broad repair/protective median AUC (pooled) | 0.940 | 0.967 |
| 4.3 cytotoxic/removal-like median AUC (pooled) | 0.755 | 0.813 |
| 4.3 any-repair cell-line delta, original AUC | 0.198 | 0.178 |
| 4.3 any-repair cell-line delta, clipped AUC | 0.193 | 0.156 |
| 4.7 broad-only repair median AUC | 0.940 | 0.970 |
| 4.7 broad-only delta | 0.194 | 0.161 |
| 4.7 broad-only MW-AUC | 0.989 | 0.985 |
| 4.7 broad-minus-glucocorticoid delta | 0.148 | 0.143 |
| 4.7 broad-minus-anti-inflammatory delta | 0.247 | 0.254 |

Confirmed unaffected (rank-based or independent of PRISM Secondary AUC, so
insensitive to the stale-annotation issue above): Table 1 compound counts;
RTR Spearman correlations (rho = 0.361 / 0.679, n = 578); TP53 interaction
(p = 0.242 / 0.276); Table 3 DepMap CRISPR re-screen; Table 4 PRISM
removal-axis re-screen; Figure 4 lineage enrichment (p = 2.53e-05); ATM
lineage-adjusted OLS (p = 0.0188).

### Fixed: S9 lineage-enrichment methodology
- S9 now uses the Bone+Soft Tissue aggregate with a one-sided Fisher exact
  test, matching Methods Section 3.13, reproducing the manuscript's
  p = 2.53e-05 / FDR = 5.82e-4 exactly.

### Added: S6 broad subclass-exclusion sensitivity
- `broad_subclass_exclusion_sensitivity.csv` added (delta-of-group-medians
  convention, matching scripts 02/03): glucocorticoid-excluded = 0.143,
  anti-inflammatory-excluded = 0.254.

### S14
- File size and MD5/SHA256 checksums recorded for
  `secondary-screen-dose-response-curve-parameters.csv`.
