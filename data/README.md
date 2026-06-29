# Data directory

Large third-party raw matrices are not redistributed in the lightweight repository bundle. For full local regeneration with `python analysis/run_all.py`, place the required public raw files under `data/` or the `data/raw/...` paths recorded in `tables/S14_input_file_provenance_manifest.csv`.

Required full-regeneration files include:

- `secondary-screen-dose-response-curve-parameters.csv`
- `secondary-screen-replicate-collapsed-logfold-change.csv` or `secondary-screen-logfold-change.csv`
- `primary-screen-replicate-collapsed-logfold-change.csv`
- `primary-screen-replicate-collapsed-treatment-info.csv`
- `OmicsSomaticMutations.csv`
- `CRISPRGeneEffect.csv`
- `Model.csv`

Optional for the auxiliary combination analysis:

- `DTP_NCI60_ALMANAC_COMBO_SCORE.xlsx` or `DTP_NCI60_ALMANAC_COMBO_SCORE.zip` under `data/raw/NCI_ALMANAC/`

Manuscript-supporting output copies are kept under `outputs/tables/` or `tables/` where applicable. Analysis scripts rebuild the corresponding outputs from raw inputs during full local regeneration. See `REGENERATION_STATUS.md` for the regeneration scope.
