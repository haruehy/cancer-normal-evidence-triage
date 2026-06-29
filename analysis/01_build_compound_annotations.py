#!/usr/bin/env python3
"""Build manuscript compound annotations from a fixed curated table.

This script intentionally does NOT re-classify PRISM compounds by keyword/regex.
The manuscript categories are based on the final curated annotation table used for
analysis and reporting. Regex helpers are used only to add descriptive subclasses
and removal-axis labels when these columns are absent; they do not change the fixed
analysis_group assignment.
"""

import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

REQUIRED_COLUMNS = [
    'drug_key', 'name', 'moa', 'target', 'analysis_group',
    'strict_repair', 'broad_repair', 'cytotoxic'
]
EXPECTED_COUNTS = {
    'strict repair/protective': 6,
    'broad repair/protective': 75,
    'cytotoxic/removal-like': 157,
    'other': 1261,
}


def _as_bool(series):
    return series.astype(str).str.strip().str.lower().isin({'true', '1', 'yes'})


def main():
    # Fixed curated annotations are the source of truth for manuscript groups.
    ann = read(CURATED_COMPOUND_ANNOTATIONS).copy()

    missing = [c for c in REQUIRED_COLUMNS if c not in ann.columns]
    if missing:
        raise ValueError(
            'curated_compound_annotations.csv is missing required columns: '
            + ', '.join(missing)
        )

    # Normalize keys and group labels without changing curated assignments.
    ann['drug_key'] = ann['drug_key'].astype(str)
    ann['analysis_group'] = ann['analysis_group'].astype(str)
    ann['strict_repair'] = _as_bool(ann['strict_repair'])
    ann['broad_repair'] = _as_bool(ann['broad_repair'])
    ann['cytotoxic'] = _as_bool(ann['cytotoxic'])

    # Fill optional descriptive columns. These are NOT used to assign analysis_group.
    if 'broad_subclass' not in ann.columns:
        ann['broad_subclass'] = ann.apply(
            lambda r: broad_subclass(r) if r.analysis_group == 'broad repair/protective' else '',
            axis=1,
        )
    if 'removal_axis' not in ann.columns:
        ann['removal_axis'] = ann.apply(removal_axis, axis=1)
    if 'annotation_source' not in ann.columns:
        ann['annotation_source'] = 'fixed curated annotation table'
    if 'curation_note' not in ann.columns:
        ann['curation_note'] = 'Group assignment imported from final manuscript curation; no regex reclassification performed.'

    # Save manuscript-facing annotation table.
    sort_cols = [c for c in ['analysis_group', 'broad_subclass', 'name', 'drug_key'] if c in ann.columns]
    save(ann.sort_values(sort_cols), 'secondary_drug_compound_annotations.csv')

    counts = ann.groupby('analysis_group').size().reindex(EXPECTED_COUNTS.keys(), fill_value=0).reset_index(name='n_compounds')
    counts['expected_n_compounds'] = counts['analysis_group'].map(EXPECTED_COUNTS)
    counts['matches_manuscript'] = counts['n_compounds'].eq(counts['expected_n_compounds'])
    save(counts, 'compound_annotation_counts.csv')

    if not counts['matches_manuscript'].all():
        raise AssertionError(
            'Curated annotation counts do not match the manuscript counts.\n'
            + counts.to_string(index=False)
        )

    print(counts.to_string(index=False))

    # Optional primary-treatment annotation. Prefer the full primary curation/proxy table
    # when present, because it carries the manuscript Primary selected-dataset
    # classification counts (28 strict, 266 broad, 215 cytotoxic). Fallback to
    # exact secondary-curation matches only when that source table is absent.
    if PRISM_PRIMARY_INFO.exists():
        pri = read(PRISM_PRIMARY_INFO)
        pri = pri.rename(columns={c: 'name' for c in pri.columns if c.lower() in ['pert_iname', 'compound_name', 'name'] and c != 'name'})
        if 'broad_id' not in pri.columns:
            pri['broad_id'] = pri.get('name', pri.index).astype(str)
        pri['drug_key'] = pri['broad_id'].astype(str)
        source_primary = ROOT/'tables'/'S1_compound_annotations_contexts_TP53_RTR_FDR'/'source_primary_drug_level_measurable_proxies.csv'
        if source_primary.exists():
            src = pd.read_csv(source_primary, low_memory=False)
            src = src.rename(columns={'eqcl_group':'analysis_group'})
            keep = ['drug_key','analysis_group','strict_repair','broad_repair','cytotoxic']
            src = src[[c for c in keep if c in src.columns]].drop_duplicates('drug_key')
            pa = pri.merge(src, on='drug_key', how='left')
        else:
            fixed_cols = ['drug_key', 'name', 'analysis_group', 'strict_repair', 'broad_repair', 'cytotoxic', 'broad_subclass', 'removal_axis']
            fixed = ann[[c for c in fixed_cols if c in ann.columns]].drop_duplicates('drug_key')
            pa = pri.merge(fixed.drop(columns=['name'], errors='ignore'), on='drug_key', how='left')
            if 'name' in pa.columns:
                fixed_name = ann[[c for c in fixed_cols if c in ann.columns]].copy()
                fixed_name['_name_key'] = fixed_name['name'].astype(str).str.lower().str.strip()
                pa['_name_key'] = pa['name'].astype(str).str.lower().str.strip()
                missing_mask = pa['analysis_group'].isna()
                fallback = pa.loc[missing_mask, ['_name_key']].merge(
                    fixed_name.drop_duplicates('_name_key').drop(columns=['drug_key'], errors='ignore'),
                    on='_name_key', how='left', suffixes=('', '_fixed')
                )
                for col in ['analysis_group', 'strict_repair', 'broad_repair', 'cytotoxic', 'broad_subclass', 'removal_axis']:
                    if col in fallback.columns:
                        pa.loc[missing_mask, col] = fallback[col].values
                pa = pa.drop(columns=['_name_key'], errors='ignore')
        pa['analysis_group'] = pa['analysis_group'].fillna('other')
        pa['strict_repair'] = pa.get('strict_repair', False).fillna(False).astype(bool)
        pa['broad_repair'] = pa.get('broad_repair', False).fillna(False).astype(bool)
        pa['cytotoxic'] = pa.get('cytotoxic', False).fillna(False).astype(bool)
        pa['broad_subclass'] = pa.apply(lambda r: broad_subclass(r) if r.analysis_group == 'broad repair/protective' else '', axis=1)
        pa['removal_axis'] = pa.apply(removal_axis, axis=1)
        save(pa, 'primary_drug_compound_annotations_selected.csv')
        pcnt = pa.groupby('analysis_group').size().reset_index(name='n_compounds')
        save(pcnt, 'primary_compound_annotation_counts.csv')


if __name__ == '__main__':
    main()
