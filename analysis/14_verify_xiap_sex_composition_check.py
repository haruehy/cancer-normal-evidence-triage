#!/usr/bin/env python3
"""Rebuild the XIAP sex-composition check from RTR, Model.csv, and CRISPRGeneEffect.csv.

This script selects the top/bottom 20% combined-RTR groups (116 each), excludes
Bone/Soft Tissue, joins DepMap sex metadata and XIAP gene effect, and writes the
per-cell-line and summary outputs used for Repository Table S15.
"""
import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
from scipy.stats import fisher_exact, chi2_contingency


def _gene_col(df, gene):
    for c in df.columns:
        if norm_gene(c) == gene:
            return c
    raise KeyError(f'{gene} not found in CRISPRGeneEffect.csv')


def main():
    rtr=read(TABLES/'cell_line_RTR_scores_with_quadrants_and_variants.csv')
    model=read(MODEL); mid='ModelID' if 'ModelID' in model.columns else model.columns[0]; model=model.rename(columns={mid:'depmap_id'})
    ge=read(CRISPR).rename(columns={pd.read_csv(CRISPR,nrows=0).columns[0]:'depmap_id'})
    xcol=_gene_col(ge, 'XIAP')
    ge=ge[['depmap_id', xcol]].rename(columns={xcol:'XIAP_gene_effect'})
    ge['XIAP_dependency_score']=-pd.to_numeric(ge['XIAP_gene_effect'], errors='coerce')
    rtr=rtr.sort_values('RTR_score_combined_rank_0_1', ascending=False).copy()
    n=int(round(len(rtr)*0.20))
    top=set(rtr.head(n).depmap_id.astype(str)); bottom=set(rtr.tail(n).depmap_id.astype(str))
    use=rtr[rtr.depmap_id.astype(str).isin(top|bottom)].copy()
    use['group']=np.where(use.depmap_id.astype(str).isin(top),'top','bottom')
    if 'OncotreeLineage' in use.columns:
        use=use[~use.OncotreeLineage.apply(is_bone_soft)].copy()
    mcols=['depmap_id']+[c for c in ['Sex','OncotreeLineage','CCLEName','ModelID'] if c in model.columns and c!='depmap_id']
    use=use.merge(model[mcols].drop_duplicates('depmap_id'), on='depmap_id', how='left', suffixes=('','_model')).merge(ge, on='depmap_id', how='left')
    if 'Sex' not in use.columns: use['Sex']='Unknown'
    use['Sex']=use['Sex'].fillna('Unknown').replace({'': 'Unknown'})
    per=use[['depmap_id','group','RTR_score_combined_rank_0_1','combined_RTR_quadrant','OncotreeLineage','Sex','XIAP_gene_effect','XIAP_dependency_score']]
    save(per, 'S15_XIAP_sex_composition_check_per_cellline.csv')
    (ROOT/'tables'/'S15_XIAP_sex_composition_check').mkdir(parents=True, exist_ok=True)
    per.to_csv(ROOT/'tables'/'S15_XIAP_sex_composition_check'/'S15_XIAP_sex_composition_check_per_cellline.csv', index=False)
    rows=[]
    for subset_name,df in [('all_included_after_lineage_exclusion',per),('XIAP_evaluable',per[per.XIAP_dependency_score.notna()])]:
        counts=df.groupby(['group','Sex']).size().unstack(fill_value=0)
        for g in ['top','bottom']:
            if g not in counts.index: counts.loc[g]=0
        female_male=[[int(counts.loc['top'].get('Female',0)), int(counts.loc['top'].get('Male',0))], [int(counts.loc['bottom'].get('Female',0)), int(counts.loc['bottom'].get('Male',0))]]
        try: orv,fp=fisher_exact(female_male, alternative='two-sided')
        except Exception: orv,fp=np.nan,np.nan
        allcats=sorted(set(df.Sex.dropna().astype(str)))
        mat=np.array([[int(counts.loc[g].get(c,0)) for c in allcats] for g in ['top','bottom']]) if allcats else np.zeros((2,0))
        try: chi2,chip,_,_=chi2_contingency(mat)
        except Exception: chi2,chip=np.nan,np.nan
        rows.append(dict(subset=subset_name, top_n=int((df.group=='top').sum()), bottom_n=int((df.group=='bottom').sum()), top_female=int(counts.loc['top'].get('Female',0)), top_male=int(counts.loc['top'].get('Male',0)), top_unknown=int(counts.loc['top'].get('Unknown',0)), bottom_female=int(counts.loc['bottom'].get('Female',0)), bottom_male=int(counts.loc['bottom'].get('Male',0)), bottom_unknown=int(counts.loc['bottom'].get('Unknown',0)), fisher_female_male_odds_ratio=orv, fisher_female_male_p=fp, chi_square_all_sex_categories_p=chip))
    d=per[per.XIAP_dependency_score.notna()]
    xt=d.loc[d.group.eq('top'),'XIAP_dependency_score']; xb=d.loc[d.group.eq('bottom'),'XIAP_dependency_score']; mw,p=mw_auc_p(xt, xb)
    rows.append(dict(subset='XIAP_dependency_all_sex', top_n=len(xt), bottom_n=len(xb), delta_top_minus_bottom=float(xt.median()-xb.median()), MW_AUC_top_greater=mw, p=p))
    known=d[d.Sex.isin(['Male','Female'])]
    reg=ols_group_lineage(known.XIAP_dependency_score, known.group.eq('top').astype(float), known.Sex)
    rows.append(dict(subset='XIAP_dependency_sex_adjusted_OLS', **reg))
    for sex in ['Male','Female']:
        s=d[d.Sex.eq(sex)]; xt=s.loc[s.group.eq('top'),'XIAP_dependency_score']; xb=s.loc[s.group.eq('bottom'),'XIAP_dependency_score']; mw,p=mw_auc_p(xt,xb)
        rows.append(dict(subset=f'XIAP_dependency_{sex.lower()}_only', top_n=len(xt), bottom_n=len(xb), delta_top_minus_bottom=float(xt.median()-xb.median()) if len(xt) and len(xb) else np.nan, MW_AUC_top_greater=mw, p=p))
    summ=pd.DataFrame(rows)
    save(summ,'S15_XIAP_sex_composition_check_summary.csv')
    summ.to_csv(ROOT/'tables'/'S15_XIAP_sex_composition_check'/'S15_XIAP_sex_composition_check_summary.csv', index=False)

if __name__=='__main__': main()
