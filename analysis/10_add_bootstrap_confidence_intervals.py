#!/usr/bin/env python3
"""Add descriptive bootstrap confidence intervals to manuscript-supporting tables.

This script adds 95% percentile bootstrap confidence intervals for delta and
MW-AUC summaries where the underlying cell-line-level vectors are bundled or
available from the expected PRISM Secondary raw input.

It intentionally does not fabricate CIs from aggregate summary rows. For tables
whose underlying response vectors are not available in the bundled repository
(e.g., raw Primary logFC candidate-axis vectors), this script leaves the existing
summary unchanged.
"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

N_BOOT = 10000
SEED = 20260627


def percentile_ci(vals, low=2.5, high=97.5):
    vals = np.asarray(vals, dtype=float)
    vals = vals[~np.isnan(vals)]
    if len(vals) == 0:
        return np.nan, np.nan
    return float(np.percentile(vals, low)), float(np.percentile(vals, high))


def fast_mw_auc(x, y):
    """Return P(x > y) + 0.5 P(x == y), equivalent to U/(n_x*n_y)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    y_sorted = np.sort(y)
    left = np.searchsorted(y_sorted, x, side='left')
    right = np.searchsorted(y_sorted, x, side='right')
    return float((left + 0.5 * (right - left)).sum() / (len(x) * len(y)))


def bootstrap_two_group_delta_mw_ci(x, y, delta_direction="x_minus_y", mw_orientation="x_greater", n_boot=N_BOOT, seed=SEED):
    x = pd.Series(x).dropna().astype(float).to_numpy()
    y = pd.Series(y).dropna().astype(float).to_numpy()
    if len(x) < 2 or len(y) < 2:
        return dict(delta_ci_low=np.nan, delta_ci_high=np.nan, MW_AUC_ci_low=np.nan, MW_AUC_ci_high=np.nan)
    rng = np.random.default_rng(seed)
    delta_vals = np.empty(n_boot, dtype=float)
    mw_vals = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        xb = rng.choice(x, size=len(x), replace=True)
        yb = rng.choice(y, size=len(y), replace=True)
        if delta_direction == "x_minus_y":
            delta_vals[b] = np.median(xb) - np.median(yb)
        elif delta_direction == "y_minus_x":
            delta_vals[b] = np.median(yb) - np.median(xb)
        else:
            raise ValueError("delta_direction must be x_minus_y or y_minus_x")
        if mw_orientation == "x_greater":
            mw_vals[b] = fast_mw_auc(xb, yb)
        elif mw_orientation == "y_greater":
            mw_vals[b] = fast_mw_auc(yb, xb)
        else:
            raise ValueError("mw_orientation must be x_greater or y_greater")
    dl, dh = percentile_ci(delta_vals)
    ml, mh = percentile_ci(mw_vals)
    return dict(delta_ci_low=dl, delta_ci_high=dh, MW_AUC_ci_low=ml, MW_AUC_ci_high=mh)


def bootstrap_one_sample_median_ci(values, n_boot=N_BOOT, seed=SEED):
    values = pd.Series(values).dropna().astype(float).to_numpy()
    if len(values) < 2:
        return dict(median_ci_low=np.nan, median_ci_high=np.nan)
    rng = np.random.default_rng(seed)
    vals = np.empty(n_boot, dtype=float)
    for b in range(n_boot):
        vals[b] = np.median(rng.choice(values, size=len(values), replace=True))
    lo, hi = percentile_ci(vals)
    return dict(median_ci_low=lo, median_ci_high=hi)


def load_secondary_with_annotations():
    ann = read(TABLES/'secondary_drug_compound_annotations.csv')
    sec = read(PRISM_SECONDARY)
    sec = sec[sec.auc.notna()].copy()
    if 'passed_str_profiling' in sec.columns:
        sec = sec[sec.passed_str_profiling.astype(str).str.lower().isin(['true','1','yes'])].copy()
    sec['drug_key'] = sec.broad_id.astype(str)
    sec = sec.merge(ann[['drug_key','analysis_group','broad_subclass']], on='drug_key', how='left')
    sec.analysis_group = sec.analysis_group.fillna('other')
    return sec


def secondary_cell_metrics(sec, val='auc', repair_groups=('strict repair/protective','broad repair/protective')):
    cell = pd.DataFrame({'depmap_id': sorted(sec.depmap_id.dropna().unique())})
    cell = cell.merge(sec[sec.analysis_group.isin(repair_groups)].groupby('depmap_id')[val].median().rename('repair_median'), on='depmap_id', how='left')
    cell = cell.merge(sec[sec.analysis_group.eq('strict repair/protective')].groupby('depmap_id')[val].median().rename('strict_median'), on='depmap_id', how='left')
    cell = cell.merge(sec[sec.analysis_group.eq('broad repair/protective')].groupby('depmap_id')[val].median().rename('broad_median'), on='depmap_id', how='left')
    cell = cell.merge(sec[sec.analysis_group.eq('cytotoxic/removal-like')].groupby('depmap_id')[val].median().rename('cytotoxic_median'), on='depmap_id', how='left')
    cell['repair_minus_cytotoxic_delta'] = cell.repair_median - cell.cytotoxic_median
    return cell


def update_repair_vs_removal_tests(sec):
    path = TABLES/'prism_secondary_repair_vs_removal_tests.csv'
    df = read(path)
    cell = pd.DataFrame({'depmap_id': sorted(sec.depmap_id.dropna().unique())})
    for col, gs in {
        'strict_repair_median_auc': ['strict repair/protective'],
        'broad_repair_median_auc': ['broad repair/protective'],
        'any_repair_median_auc': ['strict repair/protective','broad repair/protective'],
        'cytotoxic_median_auc': ['cytotoxic/removal-like'],
    }.items():
        cell = cell.merge(sec[sec.analysis_group.isin(gs)].groupby('depmap_id').auc.median().rename(col), on='depmap_id', how='left')
    for i, row in df.iterrows():
        test = str(row['test'])
        if test.startswith('strict_repair_median_auc'):
            x = cell.strict_repair_median_auc
        elif test.startswith('broad_repair_median_auc'):
            x = cell.broad_repair_median_auc
        elif test.startswith('any_repair_median_auc'):
            x = cell.any_repair_median_auc
        else:
            continue
        ci = bootstrap_two_group_delta_mw_ci(x, cell.cytotoxic_median_auc, 'x_minus_y', 'x_greater')
        df.loc[i, 'delta_ci_low'] = ci['delta_ci_low']
        df.loc[i, 'delta_ci_high'] = ci['delta_ci_high']
        df.loc[i, 'MW_AUC_ci_low'] = ci['MW_AUC_ci_low']
        df.loc[i, 'MW_AUC_ci_high'] = ci['MW_AUC_ci_high']
    save(df, 'prism_secondary_repair_vs_removal_tests.csv')


def update_auc_clip(sec):
    out_path = TABLES/'AUC_clip_sensitivity_main_results.csv'
    df = read(out_path)
    sec = sec.copy()
    sec['auc_clip1'] = sec.auc.clip(upper=1)
    cell_by_label = {
        'original_AUC': secondary_cell_metrics(sec, 'auc'),
        'clip1_AUC': secondary_cell_metrics(sec, 'auc_clip1'),
    }
    for idx, row in df.iterrows():
        if row.get('analysis') != 'cell_line_median_repair_vs_cytotoxic':
            continue
        label = row['AUC_version']
        cell = cell_by_label[label]
        dci = bootstrap_one_sample_median_ci(cell.repair_minus_cytotoxic_delta)
        mci = bootstrap_two_group_delta_mw_ci(cell.repair_median, cell.cytotoxic_median, 'x_minus_y', 'x_greater')
        df.loc[idx, 'median_value_ci_low'] = dci['median_ci_low']
        df.loc[idx, 'median_value_ci_high'] = dci['median_ci_high']
        df.loc[idx, 'MW_AUC_ci_low'] = mci['MW_AUC_ci_low']
        df.loc[idx, 'MW_AUC_ci_high'] = mci['MW_AUC_ci_high']
    save(df, 'AUC_clip_sensitivity_main_results.csv')
    # Mirror into S6 CSV directory if present.
    s6 = ROOT/'tables'/'S6_sensitivity_analyses'/'AUC_clip_sensitivity_main_results.csv'
    if s6.exists():
        df.to_csv(s6, index=False)
        print('[written]', s6)


def update_broad_subclass(sec):
    rows = []
    scenarios = [
        ('broad_all_75', None),
        ('broad_minus_glucocorticoid_49', 'glucocorticoid/corticosteroid'),
        ('broad_minus_anti_inflammatory_36', 'anti-inflammatory/COX/LOX/NFkB/TNF'),
    ]
    cyt = sec[sec.analysis_group.eq('cytotoxic/removal-like')]
    cyt_cell = cyt.groupby('depmap_id').auc.median().rename('cytotoxic_median')
    for scenario, exclude in scenarios:
        broad = sec[sec.analysis_group.eq('broad repair/protective')].copy()
        if exclude is not None:
            broad = broad[~broad.broad_subclass.eq(exclude)].copy()
        cell = pd.DataFrame({'depmap_id': sorted(sec.depmap_id.dropna().unique())})
        cell = cell.merge(broad.groupby('depmap_id').auc.median().rename('repair_median'), on='depmap_id', how='left')
        cell = cell.merge(cyt_cell, on='depmap_id', how='left')
        mw, p = mw_auc_p(cell.repair_median, cell.cytotoxic_median)
        ci = bootstrap_two_group_delta_mw_ci(cell.repair_median, cell.cytotoxic_median, 'x_minus_y', 'x_greater')
        rows.append(dict(
            scenario=scenario,
            n_cell_lines=int(cell[['repair_median','cytotoxic_median']].dropna().shape[0]),
            median_repair=float(cell.repair_median.median()),
            median_cytotoxic=float(cell.cytotoxic_median.median()),
            delta_of_medians=float(cell.repair_median.median() - cell.cytotoxic_median.median()),
            delta_of_medians_ci_low=ci['delta_ci_low'],
            delta_of_medians_ci_high=ci['delta_ci_high'],
            MW_AUC=mw,
            MW_AUC_ci_low=ci['MW_AUC_ci_low'],
            MW_AUC_ci_high=ci['MW_AUC_ci_high'],
            p_value=p,
        ))
    df = pd.DataFrame(rows)
    # Existing one-sided values are complemented with exact recomputed two-sided p-values for this output.
    out = ROOT/'tables'/'S6_sensitivity_analyses'/'broad_subclass_exclusion_sensitivity.csv'
    df.to_csv(out, index=False)
    print('[written]', out)
    save(df, 'broad_subclass_exclusion_sensitivity.csv')



def update_auc_clip_tp53_cell_line_delta():
    """Archive clipped-Secondary TP53 cell-line delta sensitivity in S6.

    This supports the manuscript statement that AUC clipping gives a Secondary
    cell-line delta p-value of approximately 0.539 for the TP53 LoF/high-impact
    versus mutation-not-called contrast.
    """
    cell_path = ROOT/'tables'/'S6_sensitivity_analyses'/'AUC_clip_sensitivity_cell_line_level.csv'
    if not cell_path.exists():
        cell_path = TABLES/'AUC_clip_sensitivity_cell_line_level.csv'
    cell = pd.read_csv(cell_path)
    rtr = read(TABLES/'cell_line_RTR_scores_with_quadrants_and_variants.csv')[['depmap_id','TP53_refined_status','any_TP53_LoF_high_impact']]
    d = cell.merge(rtr, on='depmap_id', how='left')
    rows = []
    contrasts = [('TP53 LoF/high-impact', 'TP53 mutation-not-called', 'LoF_high_impact_vs_not_called')]
    for status, ref, label in contrasts:
        x = d.loc[d.TP53_refined_status.eq(status), 'clip1_repair_minus_cytotoxic_delta'].dropna().astype(float)
        y = d.loc[d.TP53_refined_status.eq(ref), 'clip1_repair_minus_cytotoxic_delta'].dropna().astype(float)
        if len(x) < 2 or len(y) < 2:
            mw, p = np.nan, np.nan
            diff = np.nan
            med_x = np.nan
            med_y = np.nan
        else:
            res = mannwhitneyu(x, y, alternative='two-sided', method='asymptotic')
            mw = float(res.statistic / (len(x) * len(y)))
            p = float(res.pvalue)
            med_x = float(x.median())
            med_y = float(y.median())
            diff = float(med_x - med_y)
        rows.append(dict(
            analysis='AUC_clip1_secondary_cell_line_delta_TP53_interaction_sensitivity',
            dataset='secondary',
            AUC_version='clip1_AUC',
            value_metric='AUC',
            tp53_contrast=label,
            tp53_contrast_status=status,
            reference_status=ref,
            cell_line_delta_n=int(len(x) + len(y)),
            cell_line_delta_n_tp53_contrast=int(len(x)),
            cell_line_delta_n_reference=int(len(y)),
            median_delta_tp53_contrast=med_x,
            median_delta_reference=med_y,
            cell_line_delta_difference_tp53_minus_ref=diff,
            cell_line_delta_MW_AUC_tp53_greater=mw,
            cell_line_delta_p_value=p,
            test='two-sided Mann-Whitney U on clipped Secondary cell-line repair-minus-cytotoxic deltas',
            source_cell_line_table='AUC_clip_sensitivity_cell_line_level.csv',
            source_status_table='cell_line_RTR_scores_with_quadrants_and_variants.csv',
        ))
    df = pd.DataFrame(rows)
    save(df, 'AUC_clip_TP53_cell_line_delta.csv')
    s6 = ROOT/'tables'/'S6_sensitivity_analyses'/'AUC_clip_TP53_cell_line_delta.csv'
    df.to_csv(s6, index=False)
    print('[written]', s6)


def update_tp53_ci():
    rtr = read(TABLES/'cell_line_RTR_scores_with_quadrants_and_variants.csv')
    targets = [TABLES/'TP53_interaction_tests_value_level_and_cell_line_delta.csv', ROOT/'tables'/'S1_compound_annotations_contexts_TP53_RTR_FDR'/'TP53_interaction_analysis.csv']
    for p in targets:
        if not Path(p).exists():
            continue
        df = pd.read_csv(p)
        for idx, row in df.iterrows():
            dataset = str(row.get('dataset','')).lower()
            contrast = str(row.get('tp53_contrast',''))
            if dataset == 'secondary':
                metric = 'secondary_RTR_delta_any_repair_minus_removal'
            elif dataset == 'primary':
                metric = 'primary_RTR_delta_any_repair_minus_removal'
            else:
                continue
            if contrast == 'LoF_high_impact_vs_not_called':
                xstatus = 'TP53 LoF/high-impact'
            elif contrast == 'hotspot_vs_not_called':
                xstatus = 'TP53 hotspot'
            else:
                continue
            x = rtr.loc[rtr.TP53_refined_status.eq(xstatus), metric]
            y = rtr.loc[rtr.TP53_refined_status.eq('TP53 mutation-not-called'), metric]
            ci = bootstrap_two_group_delta_mw_ci(x, y, 'x_minus_y', 'x_greater')
            df.loc[idx, 'cell_line_delta_difference_ci_low'] = ci['delta_ci_low']
            df.loc[idx, 'cell_line_delta_difference_ci_high'] = ci['delta_ci_high']
            df.loc[idx, 'cell_line_delta_MW_AUC_ci_low'] = ci['MW_AUC_ci_low']
            df.loc[idx, 'cell_line_delta_MW_AUC_ci_high'] = ci['MW_AUC_ci_high']
        df.to_csv(p, index=False)
        print('[written]', p)


def main():
    sec = load_secondary_with_annotations()
    update_repair_vs_removal_tests(sec)
    update_auc_clip(sec)
    update_broad_subclass(sec)
    update_auc_clip_tp53_cell_line_delta()
    update_tp53_ci()

if __name__ == '__main__':
    main()
