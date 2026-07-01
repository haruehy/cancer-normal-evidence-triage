#!/usr/bin/env python3
"""Run the from-raw / from-curated-input analysis pipeline.

The full pipeline does not redistribute large third-party raw matrices. Place the
raw public datasets under data/ or under the data/raw/... paths recorded in
Repository Table S14, then run:

    python analysis/run_all.py

For lightweight validation of bundled manuscript-supporting outputs without raw
public datasets, run:

    python analysis/validate_repository.py
"""
from __future__ import annotations
import subprocess, sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]

def exists_any(*rels: str) -> bool:
    return any((root / r).exists() for r in rels)

required_groups = {
    'PRISM Secondary dose-response': ['data/secondary-screen-dose-response-curve-parameters.csv','data/raw/PRISM/secondary-screen-dose-response-curve-parameters.csv'],
    'PRISM Secondary logFC': ['data/secondary-screen-replicate-collapsed-logfold-change.csv','data/secondary-screen-logfold-change.csv','data/raw/PRISM/secondary-screen-replicate-collapsed-logfold-change.csv','data/raw/PRISM/secondary-screen-logfold-change.csv'],
    'PRISM Primary logFC': ['data/primary-screen-replicate-collapsed-logfold-change.csv','data/raw/PRISM/primary-screen-replicate-collapsed-logfold-change.csv'],
    'PRISM Primary treatment info': ['data/primary-screen-replicate-collapsed-treatment-info.csv','data/raw/PRISM/primary-screen-replicate-collapsed-treatment-info.csv'],
    'DepMap Model metadata': ['data/Model.csv'],
    'DepMap CRISPR gene effect': ['data/CRISPRGeneEffect.csv'],
    'DepMap OmicsSomaticMutations': ['data/OmicsSomaticMutations.csv','data/raw/DepMap/OmicsSomaticMutations.csv'],
}
missing = [name for name, rels in required_groups.items() if not exists_any(*rels)]
if missing:
    print('Full from-raw regeneration requires raw public input datasets that are not bundled.')
    print('Missing required input groups:')
    for name in missing:
        print(f'  - {name}')
    print('\nPlace the files under data/ or the data/raw/... paths recorded in tables/S14_input_file_provenance_manifest.csv.')
    print('For bundled-file validation, run: python analysis/validate_repository.py')
    raise SystemExit(2)

scripts = [
    '01_build_compound_annotations.py',
    '02_prism_repair_vs_removal_analysis.py',
    '03_auc_clip_sensitivity.py',
    '04_rtr_score_decomposition.py',
    '05_tp53_interaction_tests.py',
    '06_depmap_crispr_rescreen.py',
    '07_prism_removal_axis_rescreen.py',
    '08_almanac_combo_analysis.py',
    '09_generate_tables_and_figures.py',
    '10_add_bootstrap_confidence_intervals.py',
    '11_add_two_sided_fisher_s9.py',
    '12_curation_consistency_audit.py',
    '13_rtr_concordance_and_primary_only_enrichment.py',
    '14_verify_xiap_sex_composition_check.py',
    '15_dual_evaluable_top100_lineage_enrichment.py',
    '16_rebuild_chembl_audit_from_curated_extract.py',
    '17_rebuild_pc_descriptor.py',
    '18_count_reconciliation_tables.py',
]

for script in scripts:
    print(f'\n=== {script} ===')
    result = subprocess.run([sys.executable, str(root / 'analysis' / script)], cwd=root)
    if result.returncode:
        raise SystemExit(result.returncode)
