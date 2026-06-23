#!/usr/bin/env python3

import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
def main():
    ge=read(CRISPR).rename(columns={pd.read_csv(CRISPR,nrows=0).columns[0]:'depmap_id'}); model=read(MODEL); mid='ModelID' if 'ModelID' in model.columns else model.columns[0]; model=model.rename(columns={mid:'depmap_id'}); rtr=read(TABLES/'cell_line_RTR_scores_with_quadrants_and_variants.csv')
    cmap={norm_gene(c):c for c in ge.columns}; genes=[g for g in TARGET_GENES if g in cmap]; df=ge[['depmap_id']+[cmap[g] for g in genes]].rename(columns={cmap[g]:g for g in genes})
    for g in genes: df[f'{g}_dependency_score']=-pd.to_numeric(df[g],errors='coerce')
    meta=[c for c in ['depmap_id','CellLineName','OncotreeLineage','OncotreePrimaryDisease'] if c in model.columns]; rcols=[c for c in ['depmap_id','RTR_high_top20pct','RTR_low_bottom20pct','combined_RTR_quadrant','RTR_high_removal_high_quadrant_combined','RTR_score_combined_rank_0_1'] if c in rtr.columns]; df=df.merge(model[meta],on='depmap_id',how='left').merge(rtr[rcols],on='depmap_id',how='inner'); df['lineage']=df.get('OncotreeLineage',''); df['exclude_bone_or_soft_tissue']=df.lineage.apply(is_bone_soft); df['risk']=df.RTR_high_removal_high_quadrant_combined.astype(str).str.lower().isin(['true','1','yes']); df['high20']=df.RTR_high_top20pct.astype(str).str.lower().isin(['true','1','yes']); df['low20']=df.RTR_low_bottom20pct.astype(str).str.lower().isin(['true','1','yes'])
    comps={'repair_high_removal_high_vs_others':lambda d:np.where(d.risk,1,np.where(d.combined_RTR_quadrant.notna(),0,np.nan)),'RTR_high_top20_vs_bottom20':lambda d:np.where(d.high20,1,np.where(d.low20,0,np.nan)),'RTR_high_top20_vs_others':lambda d:np.where(d.high20,1,0)}; gene_axis={g:a for a,gs in AXIS_GENE_MAP.items() for g in gs}
    gene_rows=[]; reg_rows=[]; bone_rows=[]; axis_rows=[]; audit=[]
    for cname,func in comps.items():
        d=df.copy(); d['group']=func(d)
        for lin,sub in d.dropna(subset=['group']).groupby('lineage'): audit.append(dict(comparison=cname,lineage=lin,excluded_by_corrected_rule=is_bone_soft(lin),n_total=len(sub),n_positive=int((sub.group==1).sum()),n_negative=int((sub.group==0).sum())))
        for g in genes:
            val=f'{g}_dependency_score'; use=d.dropna(subset=['group',val]); x=use.loc[use.group.eq(1),val]; y=use.loc[use.group.eq(0),val]; mw,p=mw_auc_p(x,y); gene_rows.append(dict(comparison=cname,gene=g,axis=gene_axis.get(g,'other'),n_positive=len(x),n_negative=len(y),median_dependency_positive=x.median(),median_dependency_negative=y.median(),dependency_delta_positive_minus_negative=x.median()-y.median(),MW_AUC_positive_more_dependent=mw,p_value=p)); reg=ols_group_lineage(use[val],use.group,use.lineage); reg_rows.append(dict(comparison=cname,gene=g,axis=gene_axis.get(g,'other'),**reg)); nb=use[~use.exclude_bone_or_soft_tissue]; xb=nb.loc[nb.group.eq(1),val]; yb=nb.loc[nb.group.eq(0),val]; mwb,pb=mw_auc_p(xb,yb); bone_rows.append(dict(comparison=cname,gene=g,axis=gene_axis.get(g,'other'),n_positive_before_exclusion=len(x),n_negative_before_exclusion=len(y),n_positive_after_exclusion=len(xb),n_negative_after_exclusion=len(yb),n_positive_excluded_bone_soft_tissue=int((use.group.eq(1)&use.exclude_bone_or_soft_tissue).sum()),n_negative_excluded_bone_soft_tissue=int((use.group.eq(0)&use.exclude_bone_or_soft_tissue).sum()),dependency_delta_before=x.median()-y.median(),dependency_delta_after_exclusion=xb.median()-yb.median(),MW_AUC_after_exclusion=mwb,p_value_after_exclusion=pb,exclusion_rule="exclude OncotreeLineage in {'Bone','Soft Tissue'}"))
        for axis,gs in AXIS_GENE_MAP.items():
            gs=[g for g in gs if g in genes]
            if gs:
                cols=[f'{g}_dependency_score' for g in gs]; d[f'{axis}_dependency_median']=d[cols].median(axis=1); use=d.dropna(subset=['group',f'{axis}_dependency_median']); x=use.loc[use.group.eq(1),f'{axis}_dependency_median']; y=use.loc[use.group.eq(0),f'{axis}_dependency_median']; mw,p=mw_auc_p(x,y); axis_rows.append(dict(comparison=cname,axis=axis,genes_in_axis=';'.join(gs),n_positive=len(x),n_negative=len(y),median_dependency_positive=x.median(),median_dependency_negative=y.median(),dependency_delta_positive_minus_negative=x.median()-y.median(),MW_AUC_positive_more_dependent=mw,p_value=p))
    for rows,pcol in [(gene_rows,'p_value'),(reg_rows,'p'),(bone_rows,'p_value_after_exclusion'),(axis_rows,'p_value')]:
        t=pd.DataFrame(rows)
        for c in t.comparison.unique():
            idx=t.comparison.eq(c); q=bh_fdr(t.loc[idx,pcol])
            for i,v in zip(t.loc[idx].index,q): rows[i]['FDR_within_comparison']=v
    save(pd.DataFrame(gene_rows),'axis_gene_dependency_summary.csv'); save(pd.DataFrame(reg_rows),'lineage_adjusted_dependency_regression.csv'); save(pd.DataFrame(bone_rows),'bone_soft_tissue_excluded_dependency_corrected.csv'); save(pd.DataFrame(audit),'bone_soft_tissue_exclusion_lineage_audit.csv'); save(pd.DataFrame(axis_rows),'axis_level_dependency_summary.csv');
    if 'KIF11' in genes: save(df[[c for c in ['depmap_id','CellLineName','lineage','OncotreePrimaryDisease','RTR_score_combined_rank_0_1','combined_RTR_quadrant','RTR_high_removal_high_quadrant_combined','KIF11','KIF11_dependency_score','exclude_bone_or_soft_tissue'] if c in df.columns]],'KIF11_dependency_by_RTR_group.csv')
if __name__=='__main__': main()
