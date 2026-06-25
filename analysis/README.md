# Cancer--Normal Evidence Triage reproducible analysis pipeline

This directory contains the numbered scripts for the PRISM/RTR/DepMap/combination analyses.

## Layout

```text
analysis/
  01_build_compound_annotations.py
  02_prism_repair_vs_removal_analysis.py
  03_auc_clip_sensitivity.py
  04_rtr_score_decomposition.py
  05_tp53_interaction_tests.py
  06_depmap_crispr_rescreen.py
  07_prism_removal_axis_rescreen.py
  08_almanac_combo_analysis.py
  09_generate_tables_and_figures.py
  config.py
  README.md
```

Put raw/public files in `data/`. Default filenames are defined in `analysis/config.py` and documented in `tables/S14_input_file_provenance_manifest.csv`.

## Combined RTR and TP53 intermediates

Scripts `04_rtr_score_decomposition.py` and `05_tp53_interaction_tests.py`
read bundled archived intermediates (`data/RTR_combined_decomposition_578.csv`
and `data/TP53_interaction_recovered.csv`) so that the manuscript values can be
reproduced without re-deriving them from the very large raw PRISM/DepMap
matrices. Script 04 reproduces the combined decomposition over n=578 cell lines
(Spearman rho = 0.361 repair-preservation / 0.679 removal-vulnerability), and
script 05 reproduces the TP53 cell-line delta results (p = 0.242 Secondary,
p = 0.276 Primary). See `data/README.md` for how to rebuild these from raw data.

## Quick bundled-file validation

From the repository root, run:

```bash
python validate_repository.py
```

This does not require large raw public datasets.

## Full run order

After raw public input files are present, run:

```bash
python analysis/01_build_compound_annotations.py
python analysis/02_prism_repair_vs_removal_analysis.py
python analysis/03_auc_clip_sensitivity.py
python analysis/04_rtr_score_decomposition.py
python analysis/05_tp53_interaction_tests.py
python analysis/06_depmap_crispr_rescreen.py
python analysis/07_prism_removal_axis_rescreen.py
python analysis/08_almanac_combo_analysis.py
python analysis/09_generate_tables_and_figures.py
```

Or:

```bash
python run_all.py
```

## Key checks included

- AUC clipping sensitivity uses `AUC_clip1 = min(AUC, 1.0)`.
- Broad repair/protective subclass annotation is centralized in `config.py`.
- DepMap CRISPR uses `dependency_score = -GeneEffect`.
- Bone/Soft Tissue exclusion is corrected as separate lineages: `{"Bone", "Soft Tissue"}`. Do not use exact string `"Bone/Soft Tissue"`.

## Compound annotation source of truth

`01_build_compound_annotations.py` imports the fixed final manuscript curation table from `data/curated_compound_annotations.csv`. It does not use regular-expression rules to assign the four main manuscript groups. This prevents the strict/broad/cytotoxic/other counts from drifting away from the curated values reported in the manuscript.
