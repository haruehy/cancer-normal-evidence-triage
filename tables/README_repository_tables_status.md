# Repository tables (S1-S14)

Every repository table referenced in the manuscript (Repository Tables S1–S14) has an
explicit, S-numbered deliverable in this folder. Multi-file tables are provided both as a
folder of CSVs and as a single consolidated `.xlsx` workbook for journal upload.
This dual packaging (per-table CSV folder plus a consolidated `.xlsx`) is
intentional: the CSV folders preserve machine-readable per-file provenance,
while the `.xlsx` workbooks are the single-file artifacts intended for journal
supplementary upload. Some processed tables also appear under both
`outputs/tables/` (the copies used for the manuscript checks in
`validate_repository.py`) and here under `tables/Sxx` (the S-numbered
supplementary deliverables); this duplication is deliberate so that the
reproducibility checks and the journal supplements each remain self-contained.
See `REPOSITORY_TABLES_INDEX.csv` for the full mapping.

| S# | Deliverable |
|----|-------------|
| S1 | `S1_compound_annotations_contexts_TP53_RTR_FDR/` + `.xlsx`
| S2 | `S2_cancer_normal_evidence_table/` + `.xlsx`
| S3 | `S3_ChEMBL_normal_like_assay_audit_csv/` + `.xlsx`
| S4 | `S4_normal_protection_evidence_literature_csv/` + `.xlsx`
| S5 | `S5_representative_RTR_high_low_cell_lines.csv`
| S6 | `S6_sensitivity_analyses/` + `.xlsx`
| S7 | `S7_pairscore_sensitivity_KIF11_fallback_TN.csv`
| S8 | `S8_pairscore_sensitivity_antagonism_penalty_A.csv`
| S9 | `S9_RTR_top100_lineage_enrichment_fisher_bh.csv`
| S10 | `S10_secondary_primary_coverage_578.csv`
| S11 | `S11_depmap_crispr_rescreen/` + `.xlsx`
| S12 | `S12_prism_removal_axis_rescreen/` + `.xlsx`
| S13 | `S13_nci_almanac_combo/` + `.xlsx`
| S14 | `S14_input_file_provenance_manifest.csv`

Key bundled outputs can be checked by running python validate_repository.py from the repository root.
