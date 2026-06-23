#!/usr/bin/env python3
"""RTR score decomposition over the combined PRISM Secondary + Primary cohort.

The combined RTR score requires BOTH PRISM screens (Secondary AUC + Primary
log-fold-change). Per-cell-line repair-preservation and removal-vulnerability
ranks are computed within each screen, averaged into the combined percentile
rank, and decomposed against the combined RTR score. This yields the manuscript
cohort of n=578 cell lines (480 in both screens + 98 Primary-only) and
Spearman rho = 0.361 (repair-preservation) / 0.679 (removal-vulnerability).

The per-cell-line combined decomposition table is archived as
data/RTR_combined_decomposition_578.csv (see data/README.md for provenance).
This script reads that archived intermediate and regenerates the published
decomposition outputs. To rebuild the archived table itself from raw PRISM
data, both the Secondary curve-parameters file and the Primary
replicate-collapsed log-fold-change matrix must be present in data/.
"""
import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

ARCHIVED_COMBINED = DATA/'RTR_combined_decomposition_578.csv'

def main():
    cell = read(ARCHIVED_COMBINED)
    cell['RTR_high_removal_high_quadrant_combined'] = cell.combined_RTR_quadrant.eq('repair-high / removal-high')
    cell['RTR_high_top20pct'] = cell.RTR_score_combined_rank_0_1 >= 0.8
    cell['RTR_low_bottom20pct'] = cell.RTR_score_combined_rank_0_1 <= 0.2
    save(cell, 'cell_line_RTR_scores_with_quadrants_and_variants.csv')

    cols = ['depmap_id','RTR_score_combined_rank_0_1','combined_RTR_quadrant',
            'combined_repair_preservation_rank','combined_removal_vulnerability_rank',
            'secondary_repair_preservation_rank','secondary_removal_vulnerability_rank',
            'primary_repair_preservation_rank','primary_removal_vulnerability_rank',
            'RTR_high_top20pct','RTR_low_bottom20pct','RTR_high_removal_high_quadrant_combined']
    cols = [c for c in cols if c in cell.columns]
    save(cell[cols], 'RTR_score_decomposed_repair_preservation_removal_vulnerability.csv')

    rows = []
    for c in ['combined_repair_preservation_rank','combined_removal_vulnerability_rank',
              'secondary_repair_preservation_rank','secondary_removal_vulnerability_rank']:
        rho, p = spearman(cell[c], cell.RTR_score_combined_rank_0_1)
        rows.append(dict(component=c, n=int(cell[c].notna().sum()),
                         median=cell[c].median(), mean=cell[c].mean(),
                         corr_with_RTR_spearman=rho, corr_with_RTR_p=p))
    df = pd.DataFrame(rows); df['fdr_corr_with_RTR'] = bh_fdr(df.corr_with_RTR_p)
    save(df, 'RTR_component_summary_and_correlations.csv')

if __name__ == '__main__':
    main()
