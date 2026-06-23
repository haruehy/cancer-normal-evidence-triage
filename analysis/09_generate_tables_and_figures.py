#!/usr/bin/env python3

import os,sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import *
def load(name):
    p=TABLES/name
    return pd.read_csv(p) if p.exists() else pd.DataFrame()
def main():
    keys=['compound_annotation_counts.csv','prism_secondary_repair_vs_removal_tests.csv','AUC_clip_sensitivity_RTR_quadrant_overlap.csv','TP53_interaction_tests_value_level_and_cell_line_delta.csv','axis_gene_dependency_summary.csv','bone_soft_tissue_excluded_dependency_corrected.csv','prism_removal_axis_summary.csv','combo_repair_protection_removal_axis_summary.csv']
    save(pd.DataFrame([dict(table=k,exists=(TABLES/k).exists(),n_rows=len(load(k))) for k in keys]),'manuscript_key_results_summary.csv')
if __name__=='__main__': main()
