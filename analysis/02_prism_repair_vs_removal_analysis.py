#!/usr/bin/env python3

import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
def main():
    ann=read(TABLES/'secondary_drug_compound_annotations.csv'); sec=read(PRISM_SECONDARY); sec=sec[sec.auc.notna()].copy()
    if 'passed_str_profiling' in sec: sec=sec[sec.passed_str_profiling.astype(str).str.lower().isin(['true','1','yes'])]
    sec['drug_key']=sec.broad_id.astype(str); sec=sec.merge(ann[['drug_key','analysis_group']],on='drug_key',how='left'); sec.analysis_group=sec.analysis_group.fillna('other')
    summ=sec.groupby('analysis_group').agg(n_values=('auc','size'),n_cell_lines=('depmap_id','nunique'),n_drugs=('drug_key','nunique'),median_auc=('auc','median'),mean_auc=('auc','mean'),frac_auc_gt1=('auc',lambda x:float((x>1).mean()))).reset_index(); save(summ,'prism_secondary_repair_vs_removal_group_summary.csv')
    cell=pd.DataFrame({'depmap_id':sorted(sec.depmap_id.dropna().unique())})
    for lab,gs in {'strict_repair_median_auc':['strict repair/protective'],'broad_repair_median_auc':['broad repair/protective'],'any_repair_median_auc':['strict repair/protective','broad repair/protective'],'cytotoxic_median_auc':['cytotoxic/removal-like']}.items(): cell=cell.merge(sec[sec.analysis_group.isin(gs)].groupby('depmap_id').auc.median().rename(lab),on='depmap_id',how='left')
    cell['repair_minus_removal_delta']=cell.any_repair_median_auc-cell.cytotoxic_median_auc; cell['repair_preservation_rank']=cell.any_repair_median_auc.rank(pct=True); cell['removal_vulnerability']=1-cell.cytotoxic_median_auc; cell['removal_vulnerability_rank']=cell.removal_vulnerability.rank(pct=True); save(cell,'prism_secondary_cell_line_context.csv')
    rows=[]
    for col in ['strict_repair_median_auc','broad_repair_median_auc','any_repair_median_auc']:
        mw,p=mw_auc_p(cell[col],cell.cytotoxic_median_auc); rows.append(dict(test=f'{col}_vs_cytotoxic',median_repair=cell[col].median(),median_cytotoxic=cell.cytotoxic_median_auc.median(),delta=cell[col].median()-cell.cytotoxic_median_auc.median(),MW_AUC_repair_greater=mw,p_value=p))
    df=pd.DataFrame(rows); df['FDR']=bh_fdr(df.p_value); save(df,'prism_secondary_repair_vs_removal_tests.csv')
if __name__=='__main__': main()
