# Analysis scripts

The main pipeline is `analysis/run_all.py`. It is intended for full local regeneration after the required raw public input files have been placed under `data/` or the `data/raw/...` paths recorded in `tables/S14_input_file_provenance_manifest.csv`.

`analysis/04_rtr_score_decomposition.py` rebuilds the combined RTR cell-line table from Secondary and Primary PRISM response inputs. `analysis/05_tp53_interaction_tests.py` derives TP53 status from `OmicsSomaticMutations.csv`. `analysis/07_prism_removal_axis_rescreen.py` rebuilds Secondary AUC, Primary logFC, and Secondary logFC removal-axis screens when the corresponding raw matrices are present. `analysis/14_verify_xiap_sex_composition_check.py` rebuilds S15 from RTR outputs, `Model.csv`, and `CRISPRGeneEffect.csv`. `analysis/17_rebuild_pc_descriptor.py` rebuilds the final four-anchor `P_C` descriptor from the refined PRISM-derived survival-preservation rank; the `binary_repair_risk_rank` field is retained only as source-table provenance and is not used in the final descriptor.

The ChEMBL normal-like toxicity audit is reproducible from the archived 3,283-row extracted-and-curated ChEMBL row table onward using `analysis/16_rebuild_chembl_audit_from_curated_extract.py`. The original live ChEMBL SQL/API/web query parameters are not reconstructed.

`analysis/validate_repository.py` is a lightweight bundled-file validator. It checks the provided manuscript-supporting files without requiring third-party raw matrices.
