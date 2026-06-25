#!/usr/bin/env python3

import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
def main():
    ann=read(TABLES/'secondary_drug_compound_annotations.csv'); rtr=read(TABLES/'cell_line_RTR_scores_with_quadrants_and_variants.csv'); sec=read(PRISM_SECONDARY); sec=sec[sec.auc.notna()].copy(); sec['drug_key']=sec.broad_id.astype(str); sec=sec.merge(ann[['drug_key','analysis_group','removal_axis']],on='drug_key',how='left').merge(rtr[['depmap_id','RTR_high_removal_high_quadrant_combined','combined_RTR_quadrant']],on='depmap_id',how='inner'); sec['risk']=sec.RTR_high_removal_high_quadrant_combined.astype(str).str.lower().isin(['true','1','yes'])
    rows=[]
    for (drug,name,axis),g in sec.groupby(['drug_key','name','removal_axis'],dropna=False):
        x=g.loc[g.risk,'auc']; y=g.loc[~g.risk,'auc']
        if len(x.dropna())>=5 and len(y.dropna())>=5:
            mw,p=mw_auc_p(y,x); rows.append(dict(drug_key=drug,name=name,removal_axis=axis or 'unassigned',n_risk=len(x),n_other=len(y),median_auc_risk=x.median(),median_auc_other=y.median(),sensitivity_delta_other_minus_risk=y.median()-x.median(),MW_AUC_other_greater_than_risk=mw,p_value=p,is_definition_compound=str(g.analysis_group.iloc[0]) in {'strict repair/protective','broad repair/protective','cytotoxic/removal-like'}))
    rank=pd.DataFrame(rows).sort_values('sensitivity_delta_other_minus_risk',ascending=False); rank['FDR']=bh_fdr(rank.p_value); save(rank,'prism_removal_axis_drug_ranking.csv'); summ=rank.groupby('removal_axis').agg(n_drugs=('drug_key','nunique'),median_delta=('sensitivity_delta_other_minus_risk','median'),positive_fraction=('sensitivity_delta_other_minus_risk',lambda x:float((x>0).mean())),min_FDR=('FDR','min'),representative_drugs=('name',lambda x:'; '.join(x.dropna().astype(str).head(8)))).reset_index(); save(summ,'prism_removal_axis_summary.csv')
if __name__=='__main__': main()
