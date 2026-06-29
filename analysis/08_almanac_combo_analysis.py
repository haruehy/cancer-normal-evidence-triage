#!/usr/bin/env python3
"""Rebuild NCI-ALMANAC/combination auxiliary screen from raw combo files.

Accepted inputs under data/:
  * NCI_ALMANAC_combo_response.csv;
  * raw/NCI_ALMANAC/DTP_NCI60_ALMANAC_COMBO_SCORE.xlsx;
  * raw/NCI_ALMANAC/DTP_NCI60_ALMANAC_COMBO_SCORE.zip containing that xlsx;
  * DrugComb_combo_response.csv as a fallback.
"""
import os,sys,zipfile,tempfile; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

REPAIR={'taurine':[r'taurine'],'trolox':[r'trolox'],'benfotiamine':[r'benfotiamine'],'ergothioneine':[r'ergothioneine'],
        'amifostine':[r'amifostine'],'dexrazoxane':[r'dexrazoxane'],'celecoxib':[r'celecoxib'],
        'glucocorticoid':[r'dexamethasone|hydrocortisone|prednisolone|methylprednisolone']}

def findcol(cols,keys):
    for k in keys:
        for c in cols:
            if k.lower() in str(c).lower(): return c
    return None

def _read_combo(path):
    path=Path(path)
    if not path.exists():
        return None
    if path.suffix.lower()=='.csv':
        return pd.read_csv(path, low_memory=False)
    if path.suffix.lower() in ['.xlsx','.xls']:
        return pd.read_excel(path)
    if path.suffix.lower()=='.zip':
        with zipfile.ZipFile(path) as z:
            members=[m for m in z.namelist() if m.lower().endswith(('.xlsx','.xls','.csv')) and 'combo' in m.lower()]
            if not members:
                members=[m for m in z.namelist() if m.lower().endswith(('.xlsx','.xls','.csv'))]
            if not members: return None
            with tempfile.TemporaryDirectory() as td:
                z.extract(members[0], td)
                return _read_combo(Path(td)/members[0])
    return None

def _load_combo():
    for p in [ALMANAC, DRUGCOMB]:
        df=_read_combo(p)
        if df is not None: return p,df
    return None,None

def main():
    path,df=_load_combo()
    if df is None:
        m=pd.DataFrame([dict(status='no_combo_file_found', expected='NCI_ALMANAC_combo_response.csv, raw/NCI_ALMANAC/DTP_NCI60_ALMANAC_COMBO_SCORE.xlsx/.zip, or DrugComb_combo_response.csv', interpretation='Exact R/P candidate x removal-axis combo validation unavailable until a public combo file is added.')])
        save(m,'combo_exact_or_axis_pair_screen.csv'); save(m,'combo_repair_protection_removal_axis_summary.csv'); return
    if not any(str(c).lower() in ['drug1','drug2','drug name','drug name.1','nsc #1 b','nsc #2 b'] for c in df.columns):
        # If an Excel export has pre-header rows, re-read is difficult after read_excel; try to infer from current columns/rows.
        pass
    d1=findcol(df.columns,['drug1','drug_a','compound1','agent1','drug name','nsc #1'])
    d2=findcol(df.columns,['drug2','drug_b','compound2','agent2','drug name.1','nsc #2'])
    if d1 is None and 'Drug name' in df.columns: d1='Drug name'
    if d2 is None and 'Drug name.1' in df.columns: d2='Drug name.1'
    if not d1 or not d2:
        out=pd.DataFrame([dict(status='skipped_could_not_infer_drug_pair_columns', source_file=str(path), columns=';'.join(map(str,df.columns[:50])))])
        save(out,'combo_exact_or_axis_pair_screen.csv'); save(out,'combo_repair_protection_removal_axis_summary.csv'); return
    score=[c for c in df.columns if any(k in str(c).lower() for k in ['bliss','zip','loewe','hsa','synergy','response','viability','comboscore','combo score'])]
    meta={d1,d2,'NSC #1 b','NSC #2 b','FDA status','FDA status.1','Mechanism of action c','Mechanism of action c.1'}
    numeric_cols=[c for c in df.columns if c not in meta and pd.to_numeric(df[c],errors='coerce').notna().sum()>=10]
    if numeric_cols:
        df['median_combo_score']=df[numeric_cols].apply(pd.to_numeric,errors='coerce').median(axis=1)
        df['mean_combo_score']=df[numeric_cols].apply(pd.to_numeric,errors='coerce').mean(axis=1)
        score=list(dict.fromkeys(score+['median_combo_score','mean_combo_score']))
    rows=[]
    for _,r in df.iterrows():
        t=' | '.join([str(r.get(d1,'')),str(r.get(d2,''))])
        repair=';'.join(k for k,p in REPAIR.items() if matches(t,p))
        trap=';'.join(a for a,p in REMOVAL_AXIS_PATTERNS.items() if matches(t,p))
        if repair or trap:
            out={'drug1':r.get(d1),'drug2':r.get(d2),'repair_hit':repair,'removal_axis_hit':trap,'source_file':str(path)}
            out.update({c:r.get(c) for c in score})
            rows.append(out)
    out=pd.DataFrame(rows) if rows else pd.DataFrame([dict(status='no_exact_or_axis_pairs_found', source_file=str(path))])
    save(out,'combo_exact_or_axis_pair_screen.csv')
    if 'repair_hit' in out.columns:
        summ=out.groupby(['repair_hit','removal_axis_hit'],dropna=False).agg(n_rows=('drug1','size'))
        for c in ['median_combo_score','mean_combo_score']:
            if c in out.columns: summ[c]=out.groupby(['repair_hit','removal_axis_hit'],dropna=False)[c].median()
        summ=summ.reset_index()
    else:
        summ=out
    save(summ,'combo_repair_protection_removal_axis_summary.csv')
if __name__=='__main__': main()
