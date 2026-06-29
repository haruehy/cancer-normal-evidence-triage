# ChEMBL normal-like assay curation audit and exclusion counts

This package contains the ChEMBL normal-like assay curation/audit tables used for the manuscript.

## Main files

- `../S3_Table_ChEMBL_normal_like_assay_audit.xlsx`
  - Workbook summary with README, compound summary, accepted assays, excluded-assay sample, and curation rules.
- `chembl_normal_like_assays_curated_accepted.csv`
  - Rows accepted as normal-like/supportive assays after curation.
- `chembl_normal_like_assays_curated_excluded.csv`
  - Full excluded/ambiguous rows with `curation_reason`.
- `chembl_normal_like_assays_curated_all_rows.csv`
  - Combined accepted and excluded table used as the archived 3,283-row extraction-and-curation record.
- `chembl_curated_normal_toxicity_summary_by_compound.csv`
  - Compound-level summary used for TN prior/normal-like toxicity interpretation.
- `curation_rules.csv`
  - Curation rules.

## Exclusion-count tables

- `chembl_exclusion_reason_counts_multilabel.csv`
  - Multi-label count of each exclusion reason. A row with multiple semicolon-separated reasons contributes to each reason.
- `chembl_exclusion_reason_combination_counts.csv`
  - Counts for exact semicolon-separated reason combinations.
- `chembl_exclusion_reason_counts_by_compound.csv`
  - Multi-label exclusion reason counts by compound.

## Important interpretation note

Exclusion reasons are not mutually exclusive. Therefore, the multi-label reason counts can sum to more than the number of excluded rows.

Total excluded rows in the full excluded table: 3117
Total accepted rows in the accepted table: 166

## Extraction provenance scope

See `chembl37_extraction_query_notes.md` for the scope of the ChEMBL37 extraction provenance. The S3 audit is reproducible from the archived 3,283-row extracted set and curation rules, but the lightweight repository does not execute a live ChEMBL SQL/API extraction and does not retain the exact upstream SQL/API/web-interface query parameters.
