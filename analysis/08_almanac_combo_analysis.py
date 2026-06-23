#!/usr/bin/env python3
import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *

REPAIR={'taurine':[r'taurine'],'trolox':[r'trolox'],'benfotiamine':[r'benfotiamine'],'ergothioneine':[r'ergothioneine'],
        'amifostine':[r'amifostine'],'dexrazoxane':[r'dexrazoxane'],'celecoxib':[r'celecoxib'],
        'glucocorticoid':[r'dexamethasone|hydrocortisone|prednisolone|methylprednisolone']}

def findcol(cols,keys):
    low={str(c).lower():c for c in cols}
    for k in keys:
        for c in cols:
            if k.lower() in str(c).lower(): return c
    return None

def main():
    path=ALMANAC if ALMANAC.exists() else (DRUGCOMB if DRUGCOMB.exists() else None)
    if path is None:
        m=pd.DataFrame([dict(status='no_combo_file_found',expected='NCI_ALMANAC_combo_response.csv or DrugComb_combo_response.csv',interpretation='Exact R/P candidate x removal-axis combo validation unavailable until a public combo file is added.')])
        save(m,'combo_exact_or_axis_pair_screen.csv'); save(m,'combo_repair_protection_removal_axis_summary.csv'); return
    df=read(path)
    # Drop CellMiner metadata rows if present, or recover if the file was exported with metadata as rows.
    if not any(str(c).lower() in ['drug1','drug2','drug name','drug name.1','nsc #1 b','nsc #2 b'] for c in df.columns):
        # If the first rows contain the real header, find the row with NSC #1.
        raw=pd.read_csv(path, header=None, low_memory=False)
        header_idx=None
        for i in range(min(30,len(raw))):
            vals=[str(v) for v in raw.iloc[i].tolist()]
            if any('NSC #1' in v for v in vals) and any('NSC #2' in v for v in vals):
                header_idx=i; break
        if header_idx is not None:
            df=pd.read_csv(path, header=header_idx, low_memory=False)
    d1=findcol(df.columns,['drug1','drug_a','compound1','agent1'])
    d2=findcol(df.columns,['drug2','drug_b','compound2','agent2'])
    if d1 is None and 'Drug name' in df.columns: d1='Drug name'
    if d2 is None and 'Drug name.1' in df.columns: d2='Drug name.1'
    if not d1 or not d2:
        out=pd.DataFrame([dict(status='skipped_could_not_infer_drug_pair_columns',columns=';'.join(map(str,df.columns[:20])))])
        save(out,'combo_exact_or_axis_pair_screen.csv'); save(out,'combo_repair_protection_removal_axis_summary.csv'); return
    # score columns: named synergy/response columns, or NCI-60 numeric cell-line combo-score columns.
    score=[c for c in df.columns if any(k in str(c).lower() for k in ['bliss','zip','loewe','hsa','synergy','response','viability','comboscore','combo score'])]
    meta={d1,d2,'NSC #1 b','NSC #2 b','FDA status','FDA status.1','Mechanism of action c','Mechanism of action c.1'}
    numeric_cols=[c for c in df.columns if c not in meta and pd.to_numeric(df[c],errors='coerce').notna().sum()>=10]
    if not score and numeric_cols:
        score=numeric_cols
        df['median_combo_score']=df[numeric_cols].apply(pd.to_numeric,errors='coerce').median(axis=1)
        df['mean_combo_score']=df[numeric_cols].apply(pd.to_numeric,errors='coerce').mean(axis=1)
        score=['median_combo_score','mean_combo_score']
    rows=[]
    for _,r in df.iterrows():
        t=' | '.join([str(r.get(d1,'')),str(r.get(d2,''))])
        repair=';'.join(k for k,p in REPAIR.items() if matches(t,p))
        trap=';'.join(a for a,p in REMOVAL_AXIS_PATTERNS.items() if matches(t,p))
        if repair or trap:
            out={'drug1':r.get(d1),'drug2':r.get(d2),'repair_hit':repair,'removal_axis_hit':trap}
            out.update({c:r.get(c) for c in score})
            rows.append(out)
    out=pd.DataFrame(rows) if rows else pd.DataFrame([dict(status='no_exact_or_axis_pairs_found')])
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
