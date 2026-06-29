#!/usr/bin/env python3
"""Rebuild ChEMBL normal-like toxicity audit tables from the archived 3,283-row extract.

The live ChEMBL SQL/API query that produced the 3,283-row extract is not
reconstructed here. Starting from the archived extracted-and-curated row table,
this script deterministically regenerates accepted/excluded row tables and the
exclusion-count summaries used in Repository Table S3.
"""
import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

S3DIR = ROOT/'tables'/'S3_ChEMBL_normal_like_assay_audit_csv'


def _split_reasons(x):
    return [r for r in str(x).split(';') if r and r.lower() not in {'nan','none'}]


def main():
    require(CHEMBL_ALL_ROWS)
    all_rows = pd.read_csv(CHEMBL_ALL_ROWS, low_memory=False)
    accepted = all_rows[all_rows['curation_decision'].astype(str).str.lower().eq('accepted')].copy()
    excluded = all_rows[all_rows['curation_decision'].astype(str).str.lower().eq('excluded')].copy()
    S3DIR.mkdir(parents=True, exist_ok=True)
    accepted.to_csv(S3DIR/'chembl_normal_like_assays_curated_accepted.csv', index=False)
    excluded.to_csv(S3DIR/'chembl_normal_like_assays_curated_excluded.csv', index=False)
    print('[written]', S3DIR/'chembl_normal_like_assays_curated_accepted.csv')
    print('[written]', S3DIR/'chembl_normal_like_assays_curated_excluded.csv')
    # Multilabel exclusion counts
    rows=[]
    for reasons in excluded['curation_reason'].fillna('').map(_split_reasons):
        for r in reasons: rows.append(r)
    multilabel = pd.Series(rows).value_counts().rename_axis('curation_reason').reset_index(name='n_rows_with_reason')
    multilabel.to_csv(S3DIR/'chembl_exclusion_reason_counts_multilabel.csv', index=False)
    # Reason combinations
    comb = excluded['curation_reason'].fillna('').value_counts().rename_axis('curation_reason_combination').reset_index(name='n_rows')
    comb.to_csv(S3DIR/'chembl_exclusion_reason_combination_counts.csv', index=False)
    # By compound/reason
    cr=[]
    for _,r in excluded.iterrows():
        for reason in _split_reasons(r.get('curation_reason','')):
            cr.append((r.get('_compound'), reason))
    by = pd.DataFrame(cr, columns=['compound','curation_reason']).value_counts().rename('n_rows_with_reason').reset_index() if cr else pd.DataFrame(columns=['compound','curation_reason','n_rows_with_reason'])
    by.to_csv(S3DIR/'chembl_exclusion_reason_counts_by_compound.csv', index=False)
    # Compound summary: preserve shipped summary when present, otherwise derive a minimal deterministic summary.
    current_summary = S3DIR/'chembl_curated_normal_toxicity_summary_by_compound.csv'
    if current_summary.exists():
        summary = pd.read_csv(current_summary, low_memory=False)
    else:
        grp = all_rows.groupby('_compound', dropna=False)
        summary = grp.agg(n_raw_rows=('activity_id','size'), n_accepted_curated=('curation_decision', lambda x: int((x.astype(str).str.lower()=='accepted').sum())), n_excluded=('curation_decision', lambda x: int((x.astype(str).str.lower()=='excluded').sum()))).reset_index().rename(columns={'_compound':'compound'})
    summary.to_csv(S3DIR/'chembl_curated_normal_toxicity_summary_by_compound.csv', index=False)
    print('[ok] ChEMBL audit rebuilt from archived 3,283-row extract: accepted=', len(accepted), 'excluded=', len(excluded))

if __name__=='__main__': main()
