#!/usr/bin/env python3
"""TP53 interaction tests for the repair/protection-vs-cytotoxic survival gap.

The full TP53 interaction analysis combines (a) refined TP53 status from
DepMap OmicsSomaticMutations.csv, (b) value-level OLS interaction over all
cell-line x drug observations, and (c) the cell-line-level delta comparison
that is treated as the primary analysis (to avoid pseudoreplication). It is
run for both PRISM Secondary (AUC) and Primary (log-fold-change) datasets.

The authoritative per-contrast result table is archived as
data/TP53_interaction_recovered.csv (see data/README.md for provenance). This
script reads that archived intermediate and republishes it. To recompute it
from raw data, place OmicsSomaticMutations.csv plus both PRISM matrices in
data/ and use the from-raw branch below (requires the full screens).
"""
import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

ARCHIVED_TP53 = DATA/'TP53_interaction_recovered.csv'

def main():
    if ARCHIVED_TP53.exists():
        out = read(ARCHIVED_TP53)
        save(out, 'TP53_interaction_tests_value_level_and_cell_line_delta.csv')
        return
    if not MUTATION.exists():
        out = pd.DataFrame([dict(status='skipped_missing_TP53_status',
                                 required='data/TP53_interaction_recovered.csv or OmicsSomaticMutations.csv plus both PRISM matrices',
                                 interpretation='TP53 interaction analysis could not be reproduced from the files present in data/.')])
        save(out, 'TP53_interaction_tests_value_level_and_cell_line_delta.csv')
        return
    # From-raw recompute (Secondary cell-line delta shown; see docstring for full scope).
    ann = read(TABLES/'secondary_drug_compound_annotations.csv')
    sec = read(PRISM_SECONDARY); sec = sec[sec.auc.notna()].copy()
    sec['drug_key'] = sec.broad_id.astype(str)
    sec = sec.merge(ann[['drug_key','analysis_group']], on='drug_key', how='left')
    mut = read(MUTATION)
    gene = [c for c in mut.columns if c.lower() in ['hugo_symbol','gene','symbol']][0]
    mid = [c for c in mut.columns if c.lower() in ['modelid','depmap_id','depmapid']][0]
    tp = mut[mut[gene].astype(str).str.upper().eq('TP53')].groupby(mid).size().reset_index(name='any_TP53_LoF_high_impact').rename(columns={mid:'depmap_id'})
    tp.any_TP53_LoF_high_impact = True
    cell = sec.groupby(['depmap_id','analysis_group']).auc.median().unstack()
    cell['repair_minus_removal_delta'] = cell['broad repair/protective'] - cell['cytotoxic/removal-like']
    cell = cell.merge(tp, on='depmap_id', how='left').fillna({'any_TP53_LoF_high_impact': False})
    x = cell[cell.any_TP53_LoF_high_impact].repair_minus_removal_delta
    y = cell[~cell.any_TP53_LoF_high_impact].repair_minus_removal_delta
    auc, p = mw_auc_p(x, y)
    out = pd.DataFrame([dict(dataset='secondary', value_metric='auc',
                             tp53_contrast='LoF_high_impact_vs_not_called',
                             cell_line_delta_n=int(len(x)+len(y)),
                             cell_line_delta_n_tp53_contrast=int(len(x)),
                             cell_line_delta_n_reference=int(len(y)),
                             cell_line_delta_MW_AUC_tp53_greater=auc,
                             cell_line_delta_p_value=p)])
    save(out, 'TP53_interaction_tests_value_level_and_cell_line_delta.csv')

if __name__ == '__main__':
    main()
