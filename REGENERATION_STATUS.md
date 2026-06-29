# Regeneration status

This repository supports two modes of use: bundled-file validation and full local regeneration after the required public raw input files are downloaded separately.

## Rebuilt from raw or curated row-level inputs when source files are present

- Combined RTR 578-cell-line decomposition: `analysis/04_rtr_score_decomposition.py` rebuilds from PRISM Secondary dose-response, PRISM Primary logFC, Primary treatment information, fixed compound annotations, Model metadata, and OmicsSomaticMutations.
- TP53 interaction tests: `analysis/05_tp53_interaction_tests.py` derives TP53 status from `OmicsSomaticMutations.csv` and recomputes Secondary and Primary TP53 interaction summaries.
- PRISM removal-axis re-screen: `analysis/07_prism_removal_axis_rescreen.py` rebuilds Secondary AUC, Primary logFC, and Secondary logFC results when the corresponding PRISM matrices are present.
- XIAP sex-composition check: `analysis/14_verify_xiap_sex_composition_check.py` rebuilds the S15 summary and per-cell-line tables from RTR outputs, Model.csv, and CRISPRGeneEffect.csv.
- ChEMBL normal-like toxicity audit: `analysis/16_rebuild_chembl_audit_from_curated_extract.py` rebuilds accepted/excluded rows and exclusion-count summaries from the archived 3,283-row ChEMBL extract.
- Final P_C descriptor: `analysis/17_rebuild_pc_descriptor.py` rebuilds the final four-anchor `P_C` descriptor from the refined PRISM-derived survival-preservation rank and writes `outputs/tables/final_pc_descriptor_rebuilt.csv`. `analysis/validate_repository.py` checks agreement with Repository Table S2.
- NCI-ALMANAC auxiliary analysis: `analysis/08_almanac_combo_analysis.py` can read the raw NCI-ALMANAC xlsx or zip when placed under `data/raw/NCI_ALMANAC/`.

## Scope limitation

- The live ChEMBL extraction query itself is not reconstructed; the ChEMBL audit is reproducible from the archived 3,283-row extract onward.
- The source files retain the `binary_repair_risk_rank` field as provenance. It is not used in the final `P_C` descriptor or in the final S2 `P_C` values.

## Raw files expected for full regeneration

- `secondary-screen-dose-response-curve-parameters.csv`
- `secondary-screen-replicate-collapsed-logfold-change.csv` or `secondary-screen-logfold-change.csv`
- `primary-screen-replicate-collapsed-logfold-change.csv`
- `primary-screen-replicate-collapsed-treatment-info.csv`
- `OmicsSomaticMutations.csv`
- `CRISPRGeneEffect.csv`
- `Model.csv`
- optional: `DTP_NCI60_ALMANAC_COMBO_SCORE.xlsx` or zip for S13
