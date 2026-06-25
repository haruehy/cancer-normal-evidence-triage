#!/usr/bin/env python3

import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
def cell_metrics(sec,val):
    cell=pd.DataFrame({'depmap_id':sorted(sec.depmap_id.dropna().unique())})
    for lab,gs in {'repair_median':['strict repair/protective','broad repair/protective'],'strict_median':['strict repair/protective'],'broad_median':['broad repair/protective'],'cytotoxic_median':['cytotoxic/removal-like']}.items(): cell=cell.merge(sec[sec.analysis_group.isin(gs)].groupby('depmap_id')[val].median().rename(lab),on='depmap_id',how='left')
    cell['repair_minus_cytotoxic_delta']=cell.repair_median-cell.cytotoxic_median; cell['repair_preservation_rank']=cell.repair_median.rank(pct=True); cell['removal_vulnerability']=1-cell.cytotoxic_median; cell['removal_vulnerability_rank']=cell.removal_vulnerability.rank(pct=True); cell['secondary_RTR_rank_proxy']=cell.repair_minus_cytotoxic_delta.rank(pct=True); cell['quadrant']=[quadrant(a,b) for a,b in zip(cell.repair_preservation_rank,cell.removal_vulnerability_rank)]; return cell
def main():
    ann=read(TABLES/'secondary_drug_compound_annotations.csv')
    sec=read(PRISM_SECONDARY)
    sec=sec[sec.auc.notna()].copy()
    if 'passed_str_profiling' in sec.columns:
        sec=sec[sec.passed_str_profiling.astype(str).str.lower().isin(['true','1','yes'])].copy()
    sec['auc_clip1']=sec.auc.clip(upper=1)
    sec['drug_key']=sec.broad_id.astype(str)
    sec=sec.merge(ann[['drug_key','analysis_group']],on='drug_key',how='left')
    sec.analysis_group=sec.analysis_group.fillna('other')
    orig=cell_metrics(sec,'auc'); clip=cell_metrics(sec,'auc_clip1')
    rows=[]
    for val,label,cell in [('auc','original_AUC',orig),('auc_clip1','clip1_AUC',clip)]:
        for g,sub in sec[sec.analysis_group.isin(['strict repair/protective','broad repair/protective','cytotoxic/removal-like'])].groupby('analysis_group'):
            rows.append(dict(analysis='group_distribution',AUC_version=label,group=g,n_values=len(sub),n_cell_lines=sub.depmap_id.nunique(),n_drugs=sub.drug_key.nunique(),median_value=sub[val].median(),mean_value=sub[val].mean(),fraction_AUC_gt1_original=float((sub.auc>1).mean()),fraction_AUC_changed_by_clip=float((sub.auc!=sub.auc_clip1).mean())))
        mw,p=mw_auc_p(cell.repair_median,cell.cytotoxic_median); rows.append(dict(analysis='cell_line_median_repair_vs_cytotoxic',AUC_version=label,group='any repair/protective vs cytotoxic/removal-like',n_cell_lines=cell[['repair_median','cytotoxic_median']].dropna().shape[0],median_value=cell.repair_minus_cytotoxic_delta.median(),mean_value=cell.repair_minus_cytotoxic_delta.mean(),MW_AUC_repair_greater=mw,p_value=p))
    main=pd.DataFrame(rows); main['FDR_within_AUC_clip_main']=bh_fdr(main.get('p_value')); save(main,'AUC_clip_sensitivity_main_results.csv')
    comp=orig.add_prefix('original_').rename(columns={'original_depmap_id':'depmap_id'}).merge(clip.add_prefix('clip1_').rename(columns={'clip1_depmap_id':'depmap_id'}),on='depmap_id',how='outer'); comp['quadrant_changed']=comp.original_quadrant!=comp.clip1_quadrant; comp['delta_change_clip_minus_original']=comp.clip1_repair_minus_cytotoxic_delta-comp.original_repair_minus_cytotoxic_delta; save(comp,'AUC_clip_sensitivity_cell_line_level.csv')
    top_o=set(orig.sort_values('secondary_RTR_rank_proxy',ascending=False).head(100).depmap_id); top_c=set(clip.sort_values('secondary_RTR_rank_proxy',ascending=False).head(100).depmap_id); hh_o=set(orig.loc[orig.quadrant.eq('repair-high / removal-high'),'depmap_id']); hh_c=set(clip.loc[clip.quadrant.eq('repair-high / removal-high'),'depmap_id']); rho,p=spearman(comp.original_repair_minus_cytotoxic_delta,comp.clip1_repair_minus_cytotoxic_delta)
    save(pd.DataFrame([dict(comparison='secondary_RTR_delta_original_vs_clip1',n_cell_lines_with_delta=comp[['original_repair_minus_cytotoxic_delta','clip1_repair_minus_cytotoxic_delta']].dropna().shape[0],spearman_delta_rho=rho,spearman_delta_p=p,top100_overlap_n=len(top_o&top_c),top100_jaccard=jaccard(top_o,top_c),repair_high_removal_high_original_n=len(hh_o),repair_high_removal_high_clip1_n=len(hh_c),repair_high_removal_high_overlap_n=len(hh_o&hh_c),repair_high_removal_high_jaccard=jaccard(hh_o,hh_c),quadrant_changed_n=int(comp.quadrant_changed.sum()),quadrant_changed_fraction=float(comp.quadrant_changed.mean()))]),'AUC_clip_sensitivity_RTR_quadrant_overlap.csv')
if __name__=='__main__': main()
