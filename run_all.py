#!/usr/bin/env python3
"""Run the full analysis pipeline after raw public datasets have been added.

For the lightweight bundled-file validation that does not require raw PRISM,
DepMap, ChEMBL, or NCI-ALMANAC files, run:

    python validate_repository.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parent

required_raw_inputs = [
    root / "data" / "secondary-screen-dose-response-curve-parameters.csv",
    root / "data" / "Model.csv",
    root / "data" / "CRISPRGeneEffect.csv",
]
missing = [p for p in required_raw_inputs if not p.exists()]
if missing:
    print("Full pipeline regeneration requires raw public input datasets that are not bundled.")
    print("Missing required files:")
    for p in missing:
        print(f"  - {p.relative_to(root)}")
    print("\nDownload the raw public datasets from their original sources and place them under data/.")
    print("For the bundled-file validation, run: python validate_repository.py")
    raise SystemExit(2)

scripts = [
    "01_build_compound_annotations.py",
    "02_prism_repair_vs_removal_analysis.py",
    "03_auc_clip_sensitivity.py",
    "04_rtr_score_decomposition.py",
    "05_tp53_interaction_tests.py",
    "06_depmap_crispr_rescreen.py",
    "07_prism_removal_axis_rescreen.py",
    "08_almanac_combo_analysis.py",
    "09_generate_tables_and_figures.py",
]

for script in scripts:
    print(f"\n=== {script} ===")
    result = subprocess.run([sys.executable, str(root / "analysis" / script)], cwd=root)
    if result.returncode:
        raise SystemExit(result.returncode)
