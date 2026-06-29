# S12 Table
PRISM drug-response re-screening of removal candidate axes (secondary AUC, secondary logFC, and primary logFC), with cross-dataset summaries.


## Secondary logFC provenance

`secondary_logfc_2p5_candidate_axis_results.csv` is an archived processed output. The source raw file is recorded in S14 as `secondary-screen-replicate-collapsed-logfold-change.csv`, but this bundled repository does not regenerate the secondary-logFC S12 rows from the raw matrix. See `secondary_logfc_provenance_note.md`.

## Conditional interpretation

The FDR values in the PRISM removal-axis re-screening files are conditional within-data screening statistics because the evaluated contexts were selected from related PRISM response information. They should not be read as independent validation p-values for biological removal mechanisms.
