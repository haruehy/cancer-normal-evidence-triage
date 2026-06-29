#!/usr/bin/env python3
"""Rebuild TP53 interaction tests from raw mutation and PRISM response inputs.

This script derives refined TP53 status from OmicsSomaticMutations.csv, then
recomputes the value-level OLS interaction and cell-line-level delta comparisons
for Secondary AUC and Primary logFC. Bootstrap confidence intervals are added
later by script 10 when desired.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

TP53_LOF = 'TP53 LoF/high-impact'
TP53_HOTSPOT = 'TP53 hotspot'
TP53_MISSENSE = 'TP53 missense non-hotspot'
TP53_INFRAME = 'TP53 in-frame'
TP53_NOT_CALLED = 'TP53 mutation-not-called'


def _norm_name(name):
    return re.sub(r'[^a-z0-9]+', '', str(name).lower())


def _pick_column(df, candidates, required=True):
    lookup = {_norm_name(c): c for c in df.columns}
    for cand in candidates:
        key = _norm_name(cand)
        if key in lookup:
            return lookup[key]
    if required:
        raise KeyError(f"None of the expected columns were found: {candidates}")
    return None


def _truthy(series):
    return series.astype(str).str.strip().str.lower().isin({'true','1','yes','y','t'})


def _contains(series, pattern):
    return series.astype(str).str.contains(pattern, case=False, regex=True, na=False)


def assign_tp53_refined_status(mutation_df, universe_ids=None):
    gene_col = _pick_column(mutation_df, ['HugoSymbol','Hugo_Symbol','Hugo Symbol','gene','symbol'])
    model_col = _pick_column(mutation_df, ['ModelID','ModelId','model_id','DepMap_ID','depmap_id','DepMapID'])
    likely_col = _pick_column(mutation_df, ['LikelyLoF','Likely LoF','TranscriptLikelyLof'], required=False)
    hotspot_col = _pick_column(mutation_df, ['Hotspot','isHotspot','is_hotspot'], required=False)
    consequence_col = _pick_column(mutation_df, ['MolecularConsequence','Molecular_Consequence','Consequence'], required=False)
    variant_class_col = _pick_column(mutation_df, ['Variant_Classification','VariantClassification'], required=False)
    protein_col = _pick_column(mutation_df, ['ProteinChange','Protein_Change','HGVSp'], required=False)
    tp53 = mutation_df[mutation_df[gene_col].astype(str).str.upper().eq('TP53')].copy()
    def classify(group):
        if likely_col and _truthy(group[likely_col]).any(): return TP53_LOF
        if hotspot_col and _truthy(group[hotspot_col]).any(): return TP53_HOTSPOT
        if consequence_col and _contains(group[consequence_col], r'missense_variant').any(): return TP53_MISSENSE
        inframe = pd.Series(False, index=group.index)
        if consequence_col: inframe = inframe | _contains(group[consequence_col], r'in[_-]?frame|inframe')
        if variant_class_col: inframe = inframe | _contains(group[variant_class_col], r'in[_-]?frame|inframe')
        if protein_col: inframe = inframe | _contains(group[protein_col], r'delins|in[_-]?frame|inframe|_[A-Za-z0-9]+dup|_[A-Za-z0-9]+del|p\.[A-Za-z]\d+_[A-Za-z]\d+(?:dup|del)')
        if inframe.any(): return TP53_INFRAME
        return 'TP53 mutation-called'
    if tp53.empty:
        status = pd.DataFrame(columns=['depmap_id','TP53_refined_status'])
    else:
        status = pd.DataFrame([{'depmap_id': depmap_id, 'TP53_refined_status': classify(g)} for depmap_id,g in tp53.groupby(model_col, sort=False)])
    if universe_ids is not None:
        universe = pd.DataFrame({'depmap_id': pd.Series(list(universe_ids), dtype=str)})
        status = universe.merge(status, on='depmap_id', how='left')
        status['TP53_refined_status'] = status['TP53_refined_status'].fillna(TP53_NOT_CALLED)
    status['any_TP53_LoF_high_impact'] = status['TP53_refined_status'].eq(TP53_LOF)
    return status


def read_mutation_for_tp53():
    header = pd.read_csv(require(MUTATION), nrows=0)
    wanted = {'HugoSymbol','Hugo_Symbol','Hugo Symbol','gene','symbol','ModelID','ModelId','model_id','DepMap_ID','depmap_id','DepMapID','LikelyLoF','Likely LoF','TranscriptLikelyLof','Hotspot','isHotspot','is_hotspot','MolecularConsequence','Molecular_Consequence','Consequence','Variant_Classification','VariantClassification','ProteinChange','Protein_Change','HGVSp'}
    wanted_norm = {_norm_name(x) for x in wanted}
    usecols = [c for c in header.columns if _norm_name(c) in wanted_norm]
    gene_col = next((c for c in usecols if _norm_name(c) in {_norm_name(x) for x in ['HugoSymbol','Hugo_Symbol','Hugo Symbol','gene','symbol']}), None)
    chunks=[]
    for chunk in pd.read_csv(MUTATION, usecols=usecols, chunksize=200000, low_memory=False):
        if gene_col is None:
            chunks.append(chunk)
        else:
            sub = chunk[chunk[gene_col].astype(str).str.upper().eq('TP53')].copy()
            if not sub.empty:
                chunks.append(sub)
    return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=usecols)


def _fit_interaction(y, repair_indicator, tp53_indicator):
    d = pd.DataFrame({'y':y, 'repair':repair_indicator, 'tp53':tp53_indicator}).dropna()
    if d.empty:
        return dict(value_level_n=0, value_level_beta_group_repair_vs_cytotoxic=np.nan, value_level_beta_tp53=np.nan, value_level_beta_interaction=np.nan, value_level_se_interaction=np.nan, value_level_t_interaction=np.nan, value_level_p_interaction=np.nan)
    X = np.column_stack([np.ones(len(d)), d['repair'].astype(float).values, d['tp53'].astype(float).values, (d['repair']*d['tp53']).astype(float).values])
    yv = d['y'].astype(float).values
    beta,*_ = np.linalg.lstsq(X, yv, rcond=None)
    resid = yv - X @ beta
    dof = max(len(yv) - X.shape[1], 1)
    sig = float((resid @ resid) / dof)
    try: cov = sig * np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError: cov = sig * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t = beta[3] / se[3] if se[3] > 0 else np.nan
    p = 2 * (1 - tdist.cdf(abs(t), dof)) if not np.isnan(t) else np.nan
    return dict(value_level_n=int(len(d)), value_level_beta_group_repair_vs_cytotoxic=float(beta[1]), value_level_beta_tp53=float(beta[2]), value_level_beta_interaction=float(beta[3]), value_level_se_interaction=float(se[3]), value_level_t_interaction=float(t), value_level_p_interaction=float(p))


def _cell_delta_stats(cell, tp53_status, dataset, value_metric, contrast_label, contrast_status):
    d = cell.merge(tp53_status[['depmap_id','TP53_refined_status','any_TP53_LoF_high_impact']], on='depmap_id', how='left')
    x = d.loc[d.TP53_refined_status.eq(contrast_status), 'repair_minus_removal_delta'].dropna().astype(float)
    y = d.loc[d.TP53_refined_status.eq(TP53_NOT_CALLED), 'repair_minus_removal_delta'].dropna().astype(float)
    mw,p = mw_auc_p(x,y)
    return dict(dataset=dataset, value_metric=value_metric, tp53_contrast=contrast_label,
                direction_note='positive interaction means repair/protection-vs-cytotoxic survival gap is larger in TP53-positive contrast',
                cell_line_delta_n=int(len(x)+len(y)), cell_line_delta_n_tp53_contrast=int(len(x)), cell_line_delta_n_reference=int(len(y)),
                cell_line_delta_difference_tp53_minus_ref=float(x.median()-y.median()) if len(x) and len(y) else np.nan,
                cell_line_delta_MW_AUC_tp53_greater=mw, cell_line_delta_p_value=p)


def _secondary_value_and_cell(tp53):
    ann = read(TABLES/'secondary_drug_compound_annotations.csv')
    sec = read(PRISM_SECONDARY)
    sec = sec[sec.auc.notna()].copy()
    if 'passed_str_profiling' in sec.columns:
        sec = sec[sec.passed_str_profiling.astype(str).str.lower().isin(['true','1','yes'])].copy()
    sec['drug_key'] = sec.broad_id.astype(str)
    sec = sec.merge(ann[['drug_key','analysis_group']], on='drug_key', how='left')
    sec = sec[sec.analysis_group.isin(['strict repair/protective','broad repair/protective','cytotoxic/removal-like'])].copy()
    sec['repair_indicator'] = sec.analysis_group.isin(['strict repair/protective','broad repair/protective']).astype(int)
    sec = sec.merge(tp53[['depmap_id','TP53_refined_status','any_TP53_LoF_high_impact']], on='depmap_id', how='left')
    rows=[]
    base_cell = pd.DataFrame({'depmap_id': sorted(sec.depmap_id.dropna().unique())})
    base_cell = base_cell.merge(sec[sec.repair_indicator.eq(1)].groupby('depmap_id').auc.median().rename('repair_median'), on='depmap_id', how='left')
    base_cell = base_cell.merge(sec[sec.analysis_group.eq('cytotoxic/removal-like')].groupby('depmap_id').auc.median().rename('cytotoxic_median'), on='depmap_id', how='left')
    base_cell['repair_minus_removal_delta'] = base_cell.repair_median - base_cell.cytotoxic_median
    for label,status in [('LoF_high_impact_vs_not_called', TP53_LOF), ('hotspot_vs_not_called', 'INTERMEDIATE_TIER')]:
        if status == 'INTERMEDIATE_TIER':
            keep_status = [TP53_NOT_CALLED, TP53_HOTSPOT, TP53_MISSENSE, TP53_INFRAME, 'TP53 mutation-called']
            use = sec[sec.TP53_refined_status.isin(keep_status)].copy()
            use['tp53_indicator'] = ~use.TP53_refined_status.eq(TP53_NOT_CALLED)
            tmp_status = tp53.copy(); tmp_status['TP53_refined_status'] = np.where(tmp_status['TP53_refined_status'].isin([TP53_HOTSPOT, TP53_MISSENSE, TP53_INFRAME, 'TP53 mutation-called']), 'INTERMEDIATE_TIER', tmp_status['TP53_refined_status'])
            rec = _fit_interaction(use.auc, use.repair_indicator, use.tp53_indicator.astype(int))
            rec.update(_cell_delta_stats(base_cell, tmp_status, 'secondary', 'auc', label, 'INTERMEDIATE_TIER'))
        else:
            use = sec[sec.TP53_refined_status.isin([status, TP53_NOT_CALLED])].copy()
            use['tp53_indicator'] = use.TP53_refined_status.eq(status).astype(int)
            rec = _fit_interaction(use.auc, use.repair_indicator, use.tp53_indicator)
            rec.update(_cell_delta_stats(base_cell, tp53, 'secondary', 'auc', label, status))
        rows.append(rec)
    return rows


def _primary_value_and_cell(tp53):
    ann = read(TABLES/'primary_drug_compound_annotations_selected.csv')
    matrix = pd.read_csv(require(PRISM_PRIMARY_MATRIX), low_memory=False)
    matrix = matrix.rename(columns={matrix.columns[0]:'depmap_id'})
    cols = ann[ann.analysis_group.isin(['strict repair/protective','broad repair/protective','cytotoxic/removal-like']) & ann.column_name.isin(matrix.columns)][['column_name','analysis_group']].drop_duplicates()
    long = matrix[['depmap_id'] + cols.column_name.tolist()].melt(id_vars='depmap_id', var_name='column_name', value_name='logfold_change')
    long['logfold_change'] = pd.to_numeric(long['logfold_change'], errors='coerce')
    long = long.dropna(subset=['logfold_change']).merge(cols, on='column_name', how='left')
    long['repair_indicator'] = long.analysis_group.isin(['strict repair/protective','broad repair/protective']).astype(int)
    long = long.merge(tp53[['depmap_id','TP53_refined_status','any_TP53_LoF_high_impact']], on='depmap_id', how='left')
    rows=[]
    cell = pd.DataFrame({'depmap_id': sorted(long.depmap_id.dropna().unique())})
    cell = cell.merge(long[long.repair_indicator.eq(1)].groupby('depmap_id').logfold_change.median().rename('repair_median'), on='depmap_id', how='left')
    cell = cell.merge(long[long.analysis_group.eq('cytotoxic/removal-like')].groupby('depmap_id').logfold_change.median().rename('cytotoxic_median'), on='depmap_id', how='left')
    cell['repair_minus_removal_delta'] = cell.repair_median - cell.cytotoxic_median
    for label,status in [('LoF_high_impact_vs_not_called', TP53_LOF), ('hotspot_vs_not_called', 'INTERMEDIATE_TIER')]:
        if status == 'INTERMEDIATE_TIER':
            keep_status = [TP53_NOT_CALLED, TP53_HOTSPOT, TP53_MISSENSE, TP53_INFRAME, 'TP53 mutation-called']
            use = long[long.TP53_refined_status.isin(keep_status)].copy()
            use['tp53_indicator'] = ~use.TP53_refined_status.eq(TP53_NOT_CALLED)
            tmp_status = tp53.copy(); tmp_status['TP53_refined_status'] = np.where(tmp_status['TP53_refined_status'].isin([TP53_HOTSPOT, TP53_MISSENSE, TP53_INFRAME, 'TP53 mutation-called']), 'INTERMEDIATE_TIER', tmp_status['TP53_refined_status'])
            rec = _fit_interaction(use.logfold_change, use.repair_indicator, use.tp53_indicator.astype(int))
            rec.update(_cell_delta_stats(cell, tmp_status, 'primary', 'logfold_change', label, 'INTERMEDIATE_TIER'))
        else:
            use = long[long.TP53_refined_status.isin([status, TP53_NOT_CALLED])].copy()
            use['tp53_indicator'] = use.TP53_refined_status.eq(status).astype(int)
            rec = _fit_interaction(use.logfold_change, use.repair_indicator, use.tp53_indicator)
            rec.update(_cell_delta_stats(cell, tp53, 'primary', 'logfold_change', label, status))
        rows.append(rec)
    return rows


def main():
    require(MUTATION)
    # Universe is taken from available PRISM matrices. Mutation-not-called is assigned only within that universe.
    universe = set()
    if PRISM_SECONDARY.exists():
        try: universe.update(pd.read_csv(PRISM_SECONDARY, usecols=['depmap_id'])['depmap_id'].astype(str).unique())
        except Exception: pass
    if PRISM_PRIMARY_MATRIX.exists():
        pm = pd.read_csv(PRISM_PRIMARY_MATRIX, usecols=[0])
        universe.update(pm.iloc[:,0].astype(str).unique())
    mut = read_mutation_for_tp53()
    tp53 = assign_tp53_refined_status(mut, universe_ids=universe)
    tp53.to_csv(TABLES/'TP53_refined_status_from_raw_mutation.csv', index=False)
    print('[written]', TABLES/'TP53_refined_status_from_raw_mutation.csv')
    rows = []
    rows.extend(_secondary_value_and_cell(tp53))
    if PRISM_PRIMARY_MATRIX.exists() and (TABLES/'primary_drug_compound_annotations_selected.csv').exists():
        rows.extend(_primary_value_and_cell(tp53))
    out = pd.DataFrame(rows)
    out['value_level_fdr_interaction_within_dataset'] = np.nan
    out['cell_line_delta_fdr_within_dataset'] = np.nan
    for ds in out['dataset'].dropna().unique():
        idx=out.dataset.eq(ds)
        out.loc[idx,'value_level_fdr_interaction_within_dataset'] = bh_fdr(out.loc[idx,'value_level_p_interaction'])
        out.loc[idx,'cell_line_delta_fdr_within_dataset'] = bh_fdr(out.loc[idx,'cell_line_delta_p_value'])
    save(out, 'TP53_interaction_tests_value_level_and_cell_line_delta.csv')

if __name__ == '__main__':
    main()
