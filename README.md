# Cancer--Normal Evidence Triage

This repository contains the analysis scripts, curated annotation tables, processed outputs, and supplementary audit files supporting the manuscript:

**Cancer--Normal Evidence Triage of Repair/Protection Candidates for Repair--Remove Separation Testing**

The repository is intended to support reproducibility review. Large public raw datasets are not redistributed here and should be obtained from their original sources.

## Repository structure

```text
analysis/        Analysis scripts
data/            Curated compound annotations and input-data notes
outputs/         Processed tables used for manuscript checks
tables/          Supplementary audit/provenance tables (S1-S14)
run_all.py       Full pipeline wrapper; requires raw public datasets
validate_repository.py  Lightweight validation of bundled repository files
requirements.txt
```

Some processed tables intentionally appear in two places: under `outputs/tables/`
(the copies exercised by the manuscript checks in `validate_repository.py`) and under
`tables/Sxx` (the S-numbered repository-table deliverables).
This duplication is deliberate so that the reproducibility checks and the repository-facing table packages each remain self-contained; it is not an accidental copy.

## Key files

```text
data/curated_compound_annotations.csv
outputs/tables/
tables/S3_Table_ChEMBL_normal_like_assay_audit.xlsx
tables/S4_Table_normal_protection_evidence_literature.xlsx
tables/S6_sensitivity_analyses/RTR_weight_sensitivity_top100_overlap.csv
tables/S14_input_file_provenance_manifest.csv
```

S3 documents the ChEMBL normal-like assay audit and exclusion-reason counts. S4 documents the literature-derived normal-protection prior scoring evidence and rationale, including the anchor-candidate half-point sensitivity summary. S6 includes RTR weight-sensitivity overlap checks. S14 documents input-file provenance, including releases, expected filenames, and checksum fields.

## Quick validation

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the lightweight repository check:

```bash
python validate_repository.py
```

This validation checks the bundled curated annotation table, processed outputs, supplementary audit/provenance files, and a smoke test of the compound-annotation script. It does not require the large raw public datasets.

## Full pipeline regeneration

Full regeneration requires the public PRISM, DepMap, ChEMBL, and NCI-ALMANAC input files to be downloaded separately and placed under the expected `data/` paths described in `data/README.md`, `analysis/config.py`, and `tables/S14_input_file_provenance_manifest.csv`.

After the raw public input files have been added, run:

```bash
python run_all.py
```

Processed outputs are written to:

```text
outputs/tables/
```

If the raw public files are absent, `run_all.py` will stop with a message listing the missing files and will direct users to `validate_repository.py` for bundled-file validation.

## Public input datasets

The manuscript uses public resources including PRISM Repurposing drug-response data, DepMap CRISPR gene-effect data, DepMap model metadata, ChEMBL assay records, and NCI-ALMANAC/CellMiner combination data. These raw public datasets should be downloaded from their original sources and remain subject to the terms of those sources.

## License

Code in this repository is released under the MIT License. Curated annotation tables and manuscript-derived supplementary tables are provided for reproducibility. Original public datasets from PRISM, DepMap, ChEMBL, and NCI-ALMANAC remain subject to the terms of their respective sources.
