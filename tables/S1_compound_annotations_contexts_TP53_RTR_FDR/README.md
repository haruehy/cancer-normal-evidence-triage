# S1 Table
Compound annotations, cytotoxic-vulnerable contexts, TP53 interaction analysis, RTR decomposition (578 cell lines), and FDR details. All files are real pipeline outputs.

Additional bundled count-audit files:

- `secondary_other_count_reconciliation.csv` reconciles the curated Secondary annotation count (`other` = 1261) with the PRISM Secondary response-summary count (`other` = 1264) by documenting that three response-observed drugs without curated annotation rows are assigned to `other` in the response-summary merge.
- `prism_primary_selected_compound_annotation_counts.csv` archives the derived Primary selected-dataset class counts (Strict R/P = 28, Broad R/P = 266, cytotoxic/removal-like = 215, other = 4177) as aggregate counts only; raw PRISM Primary response data are not redistributed.
