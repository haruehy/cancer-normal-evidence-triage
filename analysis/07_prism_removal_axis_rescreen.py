#!/usr/bin/env python3
"""Rebuild conditional PRISM removal-axis re-screen outputs from PRISM response inputs.

Secondary AUC is rebuilt from the Secondary dose-response curve-parameter file.
Primary logFC is rebuilt from the Primary matrix and treatment-info annotations.
Secondary logFC is rebuilt from the Secondary logFC matrix when that raw matrix is
present. No archived Secondary-logFC result table is read as an input.
"""
import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

COMPARISON = 'RTR_high_removal_high_vs_other_no_Bone_SoftTissue'
AXIS_CANON = {
    'KIF11/kinesin/spindle': 'KIF11_kinesin',
    'microtubule/tubulin': 'microtubule_tubulin',
    'PLK/Aurora/mitotic': 'PLK_Aurora_mitotic_kinase',
    'CHK/DDR/checkpoint': 'DDR_ATR_CHK_WEE_PARP',
    'HSP/proteostasis': 'HSP90_proteostasis',
    'BCL/apoptosis': 'apoptosis_IAP_BCL_MCL',
    'topoisomerase/DNA_damage': 'topoisomerase_DNA_damage',
    'antifolate/nucleotide': 'antifolate_nucleotide',
}


def _risk_context():
    rtr = read(TABLES/'cell_line_RTR_scores_with_quadrants_and_variants.csv')
    cols = ['depmap_id','RTR_high_removal_high_quadrant_combined','combined_RTR_quadrant','OncotreeLineage']
    cols = [c for c in cols if c in rtr.columns]
    risk = rtr[cols].copy()
    if 'OncotreeLineage' in risk.columns:
        risk = risk[~risk.OncotreeLineage.apply(is_bone_soft)].copy()
    risk['risk'] = risk.RTR_high_removal_high_quadrant_combined.astype(str).str.lower().isin(['true','1','yes'])
    return risk


def _axis_summary(rows, outfile):
    df = pd.DataFrame(rows)
    if df.empty:
        save(pd.DataFrame(), outfile)
        return df
    df['FDR_within_comparison'] = bh_fdr(df['p_value'])
    save(df, outfile)
    return df


def _axis_level_from_drug(df, outfile):
    if df.empty or 'axis' not in df.columns:
        save(pd.DataFrame(), outfile); return
    rows=[]
    for axis,g in df.groupby('axis'):
        rows.append(dict(comparison=COMPARISON, axis=axis, n_drugs=int(g.drug_key.nunique()), median_delta_sensitivity=float(g.delta_sensitivity.median()), median_MW_AUC=float(g.mw_auc_sensitivity.median()), min_p_value=float(g.p_value.min()), min_FDR_within_comparison=float(g.FDR_within_comparison.min()), representative_drugs='; '.join(g.sort_values('delta_sensitivity', ascending=False).name.dropna().astype(str).head(8))))
    save(pd.DataFrame(rows), outfile)


def _secondary_auc():
    ann = read(TABLES/'secondary_drug_compound_annotations.csv')
    risk = _risk_context()
    sec = read(PRISM_SECONDARY)
    sec = sec[sec.auc.notna()].copy()
    if 'passed_str_profiling' in sec.columns:
        sec = sec[sec.passed_str_profiling.astype(str).str.lower().isin(['true','1','yes'])].copy()
    sec['drug_key'] = sec.broad_id.astype(str)
    sec = sec.merge(ann[['drug_key','name','removal_axis']], on='drug_key', how='left').merge(risk[['depmap_id','risk']], on='depmap_id', how='inner')
    sec = sec[sec.removal_axis.notna() & ~sec.removal_axis.eq('unassigned')].copy()
    rows=[]
    for (drug,name,axes),g in sec.groupby(['drug_key','name','removal_axis'], dropna=False):
        for axis_raw in str(axes).split(';'):
            axis = AXIS_CANON.get(axis_raw, axis_raw.replace('/','_').replace(' ','_'))
            x = g.loc[g.risk, 'auc']; y = g.loc[~g.risk, 'auc']
            if len(x.dropna()) >= 5 and len(y.dropna()) >= 5:
                mw,p = mw_auc_p(y,x)
                rows.append(dict(comparison=COMPARISON, dataset='Secondary', value_metric='AUC', axis=axis, drug_key=drug, name=name, n_risky=len(x.dropna()), n_other=len(y.dropna()), median_risky=x.median(), median_other=y.median(), delta_sensitivity=y.median()-x.median(), mw_auc_sensitivity=mw, p_value=p))
    return _axis_summary(rows, 'secondary_auc_candidate_axis_results.csv')


def _matrix_drug_rows(matrix_path, ann_path, value_col_name, outfile):
    ann = read(ann_path)
    risk = _risk_context()
    matrix = pd.read_csv(require(matrix_path), low_memory=False)
    matrix = matrix.rename(columns={matrix.columns[0]: 'depmap_id'})
    if 'column_name' not in ann.columns:
        # Secondary logFC may only have broad IDs as columns; map by broad_id prefix.
        ann = ann.copy(); ann['column_name'] = ann['drug_key']
    valid = []
    for _,r in ann.iterrows():
        cname = str(r.get('column_name', ''))
        dkey = str(r.get('drug_key', r.get('broad_id', '')))
        candidates = [cname, dkey]
        match = next((c for c in candidates if c in matrix.columns), None)
        if match is None:
            # Secondary logFC columns are usually replicate-collapsed treatment IDs; use prefix fallback.
            pref = dkey + '::'
            matches_cols = [c for c in matrix.columns if str(c).startswith(pref)]
            if len(matches_cols) == 1:
                match = matches_cols[0]
        if match and str(r.get('removal_axis','')) not in ['', 'nan', 'unassigned']:
            valid.append((match, dkey, r.get('name', dkey), r.get('removal_axis')))
    rows=[]
    meta = pd.DataFrame(valid, columns=['column_name','drug_key','name','removal_axis']).drop_duplicates('column_name')
    if meta.empty:
        save(pd.DataFrame([dict(status='no_matching_removal_axis_columns_found')]), outfile); return pd.DataFrame()
    long = matrix[['depmap_id'] + meta.column_name.tolist()].melt(id_vars='depmap_id', var_name='column_name', value_name=value_col_name)
    long[value_col_name] = pd.to_numeric(long[value_col_name], errors='coerce')
    long = long.dropna(subset=[value_col_name]).merge(meta, on='column_name', how='left').merge(risk[['depmap_id','risk']], on='depmap_id', how='inner')
    for (drug,name,axes),g in long.groupby(['drug_key','name','removal_axis'], dropna=False):
        for axis_raw in str(axes).split(';'):
            axis = AXIS_CANON.get(axis_raw, axis_raw.replace('/','_').replace(' ','_'))
            x = g.loc[g.risk, value_col_name]; y = g.loc[~g.risk, value_col_name]
            if len(x.dropna()) >= 5 and len(y.dropna()) >= 5:
                # For logFC, lower values mean stronger depletion, so other-minus-risk > 0 means risk more sensitive.
                mw,p = mw_auc_p(y,x)
                rows.append(dict(comparison=COMPARISON, dataset='Primary' if 'primary' in str(matrix_path).lower() else 'Secondary', value_metric='logFC', axis=axis, drug_key=drug, name=name, n_risky=len(x.dropna()), n_other=len(y.dropna()), median_risky=x.median(), median_other=y.median(), delta_sensitivity=y.median()-x.median(), mw_auc_sensitivity=mw, p_value=p))
    return _axis_summary(rows, outfile)


def main():
    sec_auc = _secondary_auc()
    primary = _matrix_drug_rows(PRISM_PRIMARY_MATRIX, TABLES/'primary_drug_compound_annotations_selected.csv', 'logfc', 'primary_logfc_candidate_axis_results.csv') if PRISM_PRIMARY_MATRIX.exists() else pd.DataFrame()
    secondary_logfc = _matrix_drug_rows(PRISM_SECONDARY_LOGFC, TABLES/'secondary_drug_compound_annotations.csv', 'logfc', 'secondary_logfc_2p5_candidate_axis_results.csv') if PRISM_SECONDARY_LOGFC.exists() else pd.DataFrame([dict(status='missing_secondary_logfc_raw_matrix', expected=str(PRISM_SECONDARY_LOGFC))])
    if isinstance(secondary_logfc, pd.DataFrame) and not secondary_logfc.empty and 'status' in secondary_logfc.columns:
        save(secondary_logfc, 'secondary_logfc_2p5_candidate_axis_results.csv')
    _axis_level_from_drug(sec_auc, 'prism_removal_axis_summary.csv')
    # Cross-dataset summaries are descriptive unions of available axis-level outputs.
    frames=[]
    for label,df in [('secondary_auc',sec_auc),('primary_logfc',primary),('secondary_logfc',secondary_logfc)]:
        if isinstance(df,pd.DataFrame) and not df.empty and 'axis' in df.columns:
            tmp=df.groupby('axis').agg(n_drugs=('drug_key','nunique'), median_delta_sensitivity=('delta_sensitivity','median'), min_FDR=('FDR_within_comparison','min')).reset_index()
            tmp['source']=label; frames.append(tmp)
    cross = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    save(cross, 'cross_dataset_axis_summary.csv')
    save(pd.concat([d for d in [sec_auc, primary, secondary_logfc] if isinstance(d,pd.DataFrame) and 'drug_key' in d.columns], ignore_index=True) if any(isinstance(d,pd.DataFrame) and 'drug_key' in d.columns for d in [sec_auc,primary,secondary_logfc]) else pd.DataFrame(), 'cross_dataset_drug_summary.csv')

if __name__=='__main__':
    main()
