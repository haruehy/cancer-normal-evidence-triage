from pathlib import Path
import re, json, math
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr, t as tdist

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'data'
OUT = ROOT/'outputs'
TABLES = OUT/'tables'
FIGS = OUT/'figures'
for d in [OUT,TABLES,FIGS]: d.mkdir(parents=True, exist_ok=True)

def first_existing(*paths):
    """Return the first existing path; if none exists, return the first candidate.

    This lets the same scripts work with either the lightweight repository
    layout (data/) or the checksum/provenance layout recorded in S14
    (data/raw/<resource>/...).
    """
    paths = [Path(p) for p in paths]
    for p in paths:
        if p.exists():
            return p
    return paths[0]

PRISM_SECONDARY = first_existing(
    DATA/'secondary-screen-dose-response-curve-parameters.csv',
    DATA/'raw/PRISM/secondary-screen-dose-response-curve-parameters.csv',
)
PRISM_SECONDARY_LOGFC = first_existing(
    DATA/'secondary-screen-replicate-collapsed-logfold-change.csv',
    DATA/'secondary-screen-logfold-change.csv',
    DATA/'raw/PRISM/secondary-screen-replicate-collapsed-logfold-change.csv',
    DATA/'raw/PRISM/secondary-screen-logfold-change.csv',
)
CURATED_COMPOUND_ANNOTATIONS = DATA/'curated_compound_annotations.csv'
PRISM_PRIMARY_MATRIX = first_existing(
    DATA/'primary-screen-replicate-collapsed-logfold-change.csv',
    DATA/'raw/PRISM/primary-screen-replicate-collapsed-logfold-change.csv',
)
PRISM_PRIMARY_INFO = first_existing(
    DATA/'primary-screen-replicate-collapsed-treatment-info.csv',
    DATA/'raw/PRISM/primary-screen-replicate-collapsed-treatment-info.csv',
)
MODEL = DATA/'Model.csv'
MUTATION = first_existing(DATA/'OmicsSomaticMutations.csv', DATA/'raw/DepMap/OmicsSomaticMutations.csv')
CRISPR = DATA/'CRISPRGeneEffect.csv'
ALMANAC = first_existing(
    DATA/'NCI_ALMANAC_combo_response.csv',
    DATA/'raw/NCI_ALMANAC/DTP_NCI60_ALMANAC_COMBO_SCORE.zip',
    DATA/'raw/NCI_ALMANAC/DTP_NCI60_ALMANAC_COMBO_SCORE.xlsx',
)
DRUGCOMB = DATA/'DrugComb_combo_response.csv'
CHEMBL_ALL_ROWS = first_existing(
    ROOT/'tables/S3_ChEMBL_normal_like_assay_audit_csv/chembl_normal_like_assays_curated_all_rows.csv',
    DATA/'chembl_normal_like_assays_curated_all_rows.csv',
)

STRICT_REPAIR_NAMES = {'taurine','trolox','benfotiamine','l-ergothioneine','ergothioneine'}
BROAD_REPAIR_PATTERNS = [
 r'glucocorticoid|corticosteroid|dexamethasone|prednisolone|hydrocortisone|methylprednisolone|diflorasone|clocortolone',
 r'cyclooxygenase|cox\b|lipoxygenase|lox\b|anti[- ]?inflammatory|nf-?kb|ikk|tnf|prostaglandin|nsaid|ibuprofen|naproxen|diclofenac|celecoxib|indomethacin',
 r'\bppar\b|peroxisome proliferator|metabolic|benfotiamine|thiamine|rage|age',
 r'antioxidant|radical scavenger|trolox|taurine|ergothioneine|glutathione|n-acetyl|ascorbic|vitamin e',
 r'nitric oxide|nos inhibitor|no scavenger|nitric oxide synthase']
CYTOTOXIC_PATTERNS = [
 r'topoisomerase|alkylating|DNA damaging|DNA damage|PARP|platinum|anthracycline|antimetabolite',
 r'microtubule|tubulin|taxane|vinca|vincristine|vinblastine|epothilone|mitotic|spindle|kinesin|KIF11',
 r'proteasome|bortezomib|carfilzomib', r'apoptosis|BCL|MCL1|navitoclax|venetoclax',
 r'CHK|checkpoint|CHEK|ATR|ATM|WEE1', r'PLK|Aurora|AURK', r'HSP|HSP90']
BROAD_SUBCLASS_RULES = [
 ('glucocorticoid/corticosteroid', r'glucocorticoid|corticosteroid|dexamethasone|prednisolone|hydrocortisone|methylprednisolone|diflorasone|clocortolone'),
 ('anti-inflammatory/COX/LOX/NFkB/TNF', r'cyclooxygenase|cox\b|lipoxygenase|lox\b|anti[- ]?inflammatory|nf-?kb|ikk|tnf|prostaglandin|nsaid|ibuprofen|naproxen|diclofenac|celecoxib|indomethacin'),
 ('PPAR/metabolic-stress', r'\bppar\b|peroxisome proliferator|metabolic|benfotiamine|thiamine|rage|age'),
 ('antioxidant/radical-scavenger', r'antioxidant|radical scavenger|trolox|taurine|ergothioneine|glutathione|n-acetyl|ascorbic|vitamin e'),
 ('NO/NOS-related', r'nitric oxide|nos inhibitor|no scavenger|nitric oxide synthase')]
REMOVAL_AXIS_PATTERNS = {
 'KIF11/kinesin/spindle':[r'KIF11|kinesin|ispinesib|filanesib|litronesib|SB[- ]?743921'],
 'microtubule/tubulin':[r'microtubule|tubulin|vincristine|vinblastine|vindesine|vinflunine|epothilone|paclitaxel|docetaxel'],
 'PLK/Aurora/mitotic':[r'\bPLK\b|PLK1|Aurora|AURK|volasertib|alisertib|BI[- ]?2536|hesperadin'],
 'CHK/DDR/checkpoint':[r'\bCHK\b|CHEK|prexasertib|LY2606368|PF[- ]?477736|CHIR[- ]?124|ATR|ATM|WEE1'],
 'HSP/proteostasis':[r'\bHSP\b|HSP90|AUY922|alvespimycin|NMS[- ]?E973|geldanamycin'],
 'BCL/apoptosis':[r'\bBCL\b|MCL1|navitoclax|venetoclax'],
 'topoisomerase/DNA_damage':[r'topoisomerase|TOP1|TOP2|irinotecan|etoposide|doxorubicin'],
 'antifolate/nucleotide':[r'DHFR|TYMS|antifolate|pralatrexate|methotrexate|pemetrexed']}
TARGET_GENES = ['KIF11','PLK1','AURKA','AURKB','CHEK1','CHEK2','HSP90AA1','HSP90AB1','BCL2L1','MCL1','BCL2','TUBA1B','TUBB','KIFC1','ATM','XIAP']
AXIS_GENE_MAP = {'KIF11_spindle':['KIF11'],'PLK_Aurora_mitotic':['PLK1','AURKA','AURKB'],'CHK_DDR_checkpoint':['CHEK1','CHEK2','ATM'],'HSP_proteostasis':['HSP90AA1','HSP90AB1'],'BCL_apoptosis':['BCL2L1','MCL1','BCL2','XIAP'],'microtubule_spindle_support':['TUBA1B','TUBB','KIFC1']}

def require(path):
    path=Path(path)
    if not path.exists(): raise FileNotFoundError(f'Missing required file: {path}')
    return path

def read(path): return pd.read_csv(require(path), low_memory=False)
def save(df, name):
    path = TABLES/name; path.parent.mkdir(parents=True, exist_ok=True); df.to_csv(path,index=False); print('[written]',path); return path

def norm_gene(c): return re.sub(r'\s*\(\d+\)\s*$','',str(c)).strip()
def norm_lineage(x): return re.sub(r'\s+',' ',str(x or '').strip()).lower()
def is_bone_soft(lineage): return norm_lineage(lineage) in {'bone','soft tissue'}
def boolish(x): return str(x).strip().lower() in {'true','1','yes'}
def drug_text(row): return ' | '.join(str(row.get(c,'')) for c in ['name','Name','moa','MOA','target','Target','broad_id','column_name'])
def matches(text, pats): return any(re.search(p, str(text), re.I) for p in pats)
def analysis_group(row):
    txt=drug_text(row); name=str(row.get('name',row.get('Name',''))).strip().lower()
    if name in STRICT_REPAIR_NAMES: return 'strict repair/protective'
    if matches(txt, BROAD_REPAIR_PATTERNS): return 'broad repair/protective'
    if matches(txt, CYTOTOXIC_PATTERNS): return 'cytotoxic/removal-like'
    return 'other'
def broad_subclass(row):
    txt=drug_text(row)
    for lab,pat in BROAD_SUBCLASS_RULES:
        if re.search(pat,txt,re.I): return lab
    return 'other broad repair/protective'
def removal_axis(row):
    txt=drug_text(row); hits=[]
    for axis,pats in REMOVAL_AXIS_PATTERNS.items():
        if matches(txt,pats): hits.append(axis)
    return ';'.join(hits) if hits else 'unassigned'
def bh_fdr(pvals):
    p=np.asarray(pd.to_numeric(pd.Series(pvals),errors='coerce'),dtype=float); out=np.full(len(p),np.nan); m=~np.isnan(p); pm=p[m]
    if len(pm)==0: return out
    order=np.argsort(pm); ranked=pm[order]; q=ranked*len(pm)/np.arange(1,len(pm)+1); q=np.minimum.accumulate(q[::-1])[::-1]
    tmp=np.empty(len(pm)); tmp[order]=np.clip(q,0,1); out[m]=tmp; return out
def mw_auc_p(x,y):
    x=pd.Series(x).dropna().astype(float); y=pd.Series(y).dropna().astype(float)
    if len(x)<2 or len(y)<2: return np.nan,np.nan
    r=mannwhitneyu(x,y,alternative='two-sided',method='asymptotic'); return float(r.statistic/(len(x)*len(y))),float(r.pvalue)
def spearman(x,y):
    d=pd.DataFrame({'x':x,'y':y}).dropna()
    if len(d)<3: return np.nan,np.nan
    r=spearmanr(d.x,d.y); return float(r.statistic),float(r.pvalue)
def quadrant(repair_rank, removal_rank):
    if pd.isna(repair_rank) or pd.isna(removal_rank): return 'unclassified'
    if repair_rank>=0.5 and removal_rank>=0.5: return 'repair-high / removal-high'
    if repair_rank>=0.5 and removal_rank<0.5: return 'repair-high / removal-low'
    if repair_rank<0.5 and removal_rank>=0.5: return 'repair-low / removal-high'
    return 'repair-low / removal-low'
def jaccard(a,b):
    a,b=set(a),set(b); return len(a&b)/len(a|b) if a|b else np.nan
def ols_group_lineage(y, group, lineage):
    d=pd.DataFrame({'y':y,'g':group,'lin':lineage}).dropna()
    if len(d)<20 or d.g.nunique()<2: return dict(n=len(d),n_lineages=np.nan,beta_group=np.nan,se=np.nan,t=np.nan,p=np.nan,beta_group_z=np.nan)
    lins=sorted(d.lin.astype(str).unique()); X=[np.ones(len(d)), d.g.astype(float).values]
    for lin in lins[1:]: X.append((d.lin.astype(str).values==lin).astype(float))
    X=np.column_stack(X); yv=d.y.astype(float).values
    def fit(yv):
        beta,*_=np.linalg.lstsq(X,yv,rcond=None); resid=yv-X@beta; dof=max(len(yv)-X.shape[1],1); sig=float((resid@resid)/dof)
        try: cov=sig*np.linalg.inv(X.T@X)
        except np.linalg.LinAlgError: cov=sig*np.linalg.pinv(X.T@X)
        se=np.sqrt(np.diag(cov)); t=beta[1]/se[1] if se[1]>0 else np.nan; p=2*(1-tdist.cdf(abs(t),dof)) if not np.isnan(t) else np.nan
        return beta[1],se[1],t,p
    b,se,t,p=fit(yv); yz=(yv-yv.mean())/yv.std(ddof=0) if yv.std(ddof=0)>0 else yv*np.nan; bz=fit(yz)[0] if not np.isnan(yz).all() else np.nan
    return dict(n=len(d),n_lineages=len(lins),beta_group=b,se=se,t=t,p=p,beta_group_z=bz)
