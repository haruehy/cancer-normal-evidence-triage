# ChEMBL normal-like assay curation audit and exclusion counts

This package contains the ChEMBL normal-like assay curation/audit tables used for the manuscript.

## Main files

- `../S3_Table_ChEMBL_normal_like_assay_audit.xlsx`
  - Workbook summary with README, compound summary, accepted assays, excluded-assay sample, and curation rules.
- `chembl_normal_like_assays_curated_accepted.csv`
  - Rows accepted as normal-like/supportive assays after curation.
- `chembl_normal_like_assays_curated_excluded.csv`
  - Full excluded/ambiguous rows with `curation_reason`.
- `chembl_curated_normal_toxicity_summary_by_compound.csv`
  - Compound-level summary used for TN prior/normal-like toxicity interpretation.
- `curation_rules.csv`
  - Curation rules.

## Added exclusion-count tables

- `chembl_exclusion_reason_counts_multilabel.csv`
  - Multi-label count of each exclusion reason. A row with multiple semicolon-separated reasons contributes to each reason.
- `chembl_exclusion_reason_combination_counts.csv`
  - Counts for exact semicolon-separated reason combinations.
- `chembl_exclusion_reason_counts_by_compound.csv`
  - Multi-label exclusion reason counts by compound.

## Note on the combined (all-rows) table

The previously bundled `chembl_normal_like_assays_curated_all_rows.csv` was the
simple concatenation of the accepted and excluded tables (3283 rows = 166
accepted + 3117 excluded, identical 63-column schema). It has been removed to
avoid redundant storage and can be reconstructed at any time by concatenating
`chembl_normal_like_assays_curated_accepted.csv` and
`chembl_normal_like_assays_curated_excluded.csv`.

## Important interpretation note

Exclusion reasons are not mutually exclusive. Therefore, the multi-label reason counts can sum to more than the number of excluded rows.

Total excluded rows in the full excluded table: 3117
Total accepted rows in the accepted table: 166
