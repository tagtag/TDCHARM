#!/usr/bin/env python3
"""
build_supplementary_tables.py

Assembles Supplementary Tables S1-S7 of the revised manuscript into one
Excel workbook (one sheet per table, plus a README sheet), from the
tab-separated outputs that scripts 19-25 write to the working directory.
The values are copied unchanged; every table is labelled with the script
and file it comes from. Large outputs (per-gene scores, complete GO tables,
factor correlation matrices, ABC-like scores of all pairs) are not
included; they are written by the same scripts to the working directory.

The only value not read from a file is the descriptive sign test of the
original analysis (12 of 12 components; one-sided binomial P = 0.5^12).

Usage:
    python3 tools/build_supplementary_tables.py <working directory> <output.xlsx>

Requires: pandas, openpyxl
"""
import sys, os, pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter
if len(sys.argv) != 3:
    sys.exit('usage: build_supplementary_tables.py <working directory> <output.xlsx>')
D=sys.argv[1]; OUT=sys.argv[2]
P='GSE303006_M23_'
def rd(f): return pd.read_csv(os.path.join(D,P+f),sep='\t')
F=Font(name='Arial',size=10); B=Font(name='Arial',size=10,bold=True); T=Font(name='Arial',size=12,bold=True); I=Font(name='Arial',size=10,italic=True)
H=PatternFill('solid',start_color='DDE4EE')
wb=Workbook(); wb.remove(wb.active)
def sheet(name,title,desc,blocks):
    ws=wb.create_sheet(name); r=1
    ws.cell(r,1,title).font=T; r+=1
    for line in desc: ws.cell(r,1,line).font=I; r+=1
    r+=1; widths={}
    for btitle,src,df in blocks:
        ws.cell(r,1,btitle).font=B; ws.cell(r,1).alignment=Alignment(wrap_text=False); r+=1
        ws.cell(r,1,'Source: '+src).font=I; r+=1
        for j,c in enumerate(df.columns,1):
            x=ws.cell(r,j,str(c)); x.font=B; x.fill=H; widths[j]=max(widths.get(j,8),min(len(str(c))+2,40))
        r+=1
        for row in df.itertuples(index=False):
            for j,v in enumerate(row,1):
                if pd.isna(v): v=None
                elif hasattr(v,'item'): v=v.item()
                x=ws.cell(r,j,v); x.font=F
                if isinstance(v,float):
                    x.number_format='0.000E+00' if (v!=0 and (abs(v)<1e-3 or abs(v)>=1e5)) else '0.0000'
                widths[j]=max(widths.get(j,8),min(len(str(v))+2 if v is not None else 4,40))
            r+=1
        r+=1
    for j,w in widths.items(): ws.column_dimensions[get_column_letter(j)].width=w
    ws.freeze_panes=None
    return ws
src=lambda s,f: f'Script {s}; {P}{f}'
S=[]
S.append(('S1','Supplementary Table S1. Stage-1 and stage-2 variance retention and reconstruction errors',
 ['Stage-1 retained fraction per modality (rank 20), cumulative retention for k = 1-20, stage-2 core energy, and retention per modality.'],
 [(n,src(19,f),rd(f)) for n,f in [('Stage 1: retention per modality','rev_stage1_energy_by_modality.tsv'),('Stage 1: cumulative retention, k = 1-20','rev_stage1_cumulative_energy.tsv'),('Stage 2: conceptual tensor and core','rev_stage2_energy.tsv'),('Stage 1, stage 2 and overall retention per modality','rev_stage2_energy_by_modality.tsv')]]))
S.append(('S2','Supplementary Table S2. Rotations within near-degenerate stage-1 subspaces',
 ['Blocks of stage-1 components with relative singular-value gap < 0.05; results of each rotation (20 Haar rotations and 45-degree rotations) and summary.'],
 [(n,src(20,f),rd(f)) for n,f in [('Near-degenerate blocks','rev20_near_degenerate_blocks.tsv'),('Rotation results','rev20_rotation_results.tsv'),('Summary','rev20_rotation_summary.tsv')]]))
S.append(('S3','Supplementary Table S3. Independent reanalysis of the two-replicate subsets',
 ['Scripts 05-08 rerun on R1+R2, R1+R3 and R2+R3. Factor matching (21b), core element at the position of core 12 (21c), and the automatic selection rule (21).',
  'Full |r| matrices, gene scores and complete GO tables are written by the same scripts to the working directory and are not included here.'],
 [(n,src(s,f),rd(f)) for n,s,f in [
  ('Best-matching subset factor for every full-data factor',"21b",'rev21b_best_matches.tsv'),
  ('Core element at the position of core 12 in each run',"21c",'rev21c_core_by_position.tsv'),
  ('Gene-ranking agreement (core-12 position)',"21c",'rev21c_gene_ranking_agreement.tsv'),
  ('GO-BP agreement (core-12 position)',"21c",'rev21c_GO_agreement.tsv'),
  ('Top 10 GO-BP terms per run and side (core-12 position)',"21c",'rev21c_GO_top10_per_run.tsv'),
  ('Element selected by the automatic rule in each run',"21",'rev21_selection.tsv'),
  ('Candidate elements of the automatic rule',"21",'rev21_selection_core_candidates.tsv'),
  ('Gene-ranking agreement (rule-selected element)',"21",'rev21_gene_ranking_agreement.tsv'),
  ('GO-BP agreement (rule-selected element)',"21",'rev21_GO_agreement.tsv'),
  ('Cell-side statistics per run',"21",'rev21_cell_side.tsv')]]))
sign=pd.DataFrame({'quantity':['components among top-100 core elements','components with celltype epsilon2 > replicate epsilon2','one-sided sign test P (binomial, p = 0.5)'],'value':[12,12,0.5**12]})
S.append(('S4','Supplementary Table S4. Cell type versus replicate',
 ['Cell-type-by-replicate composition, joint linear models (type-II partial eta2), within-cell-type replicate tests (22), the core-12 contrast within replicates (22b), and the descriptive sign test of the original analysis.',
  'Per-cell covariates are written by Script 22 to the working directory and are not included here.'],
 [(n,src(s,f),rd(f)) for n,s,f in [
  ('Composition: counts',"22",'rev22_celltype_by_replicate_counts.tsv'),
  ('Composition: row proportions',"22",'rev22_celltype_by_replicate_proportions.tsv'),
  ('Joint model: summary per component',"22",'rev22_joint_model_summary.tsv'),
  ('Joint model: partial eta2 for every term',"22",'rev22_joint_model_partial_eta2.tsv'),
  ('Within-cell-type replicate Kruskal-Wallis: summary',"22",'rev22_within_celltype_replicate_summary.tsv'),
  ('Within-cell-type replicate Kruskal-Wallis: all tests',"22",'rev22_within_celltype_replicate_KW.tsv'),
  ('CELL03 by cell type and replicate',"22b",'rev22b_summary_by_type_replicate.tsv'),
  ('Ndnf/Lamp5 vs L3/4 IT contrast within each replicate',"22b",'rev22b_contrast_within_replicate.tsv'),
  ('Pairwise replicate comparisons within cell type',"22b",'rev22b_replicate_pairs_within_type.tsv')]]
 +[('Descriptive sign test (moved from the main text)','Script 10 (original analysis): 12 of 12 components; one-sided binomial P = 0.5^12, as reported in the original manuscript',sign)]))
S.append(('S5','Supplementary Table S5. Comparison with supported enhancer-gene links',
 ['Activity-by-contact-like links from the pseudobulk data (23) and cell-type-resolved links of the original CHARM study (23b; Chen et al., Nature 2026, Supplementary Table 7).',
  'ABC-like scores of all pairs and the parsed CHARM links are written by Scripts 23 and 23b to the working directory and are not included here; the CHARM workbook itself is not redistributed.'],
 [(n,src(s,f),rd(f)) for n,s,f in [
  ('ABC-like reference sets',"23",'rev23_reference_sets.tsv'),
  ('ABC-like: top-N enrichment, CMH stratified by distance decile',"23",'rev23_CMH_topN.tsv'),
  ('ABC-like: AUPRC vs distance-matched permutation null',"23",'rev23_AUPRC.tsv'),
  ('CHARM links: coverage by the candidate pairs',"23b",'rev23b_CHARM_coverage.tsv'),
  ('CHARM reference sets',"23b",'rev23b_reference_sets.tsv'),
  ('CHARM: top-N enrichment, CMH stratified by distance decile',"23b",'rev23b_CMH_topN.tsv'),
  ('CHARM: AUPRC vs distance-matched permutation null',"23b",'rev23b_AUPRC.tsv'),
  ('CHARM links specific to Ndnf/Lamp5 vs L3/4 IT: test',"23b",'rev23b_groupA_vs_groupB_links.tsv'),
  ('CHARM-linked pairs by cell-type class: summary',"23b",'rev23b_groupA_vs_groupB_links_summary.tsv')]]))
S.append(('S6','Supplementary Table S6. Modality ablation',
 ['Stage 2 recomputed for each single modality and each leave-one-modality-out set. 24: element selected by the automatic rule; 24b: factors matched directly to EP03 and CELL03 of core 12.',
  'Full |r| matrices, gene scores and GO tables are written by the same scripts to the working directory and are not included here.'],
 [(n,src(s,f),rd(f)) for n,s,f in [
  ('Core element from the factors matching core 12',"24b",'rev24b_core12_matched.tsv'),
  ('Gene-ranking agreement with core 12 (matched factors)',"24b",'rev24b_gene_ranking_agreement.tsv'),
  ('GO-BP agreement with core 12 (matched factors)',"24b",'rev24b_GO_agreement.tsv'),
  ('Best-matching factor for every four-modality factor',"24b",'rev24b_best_matches.tsv'),
  ('Variants: element selected by the automatic rule',"24",'rev24_variants.tsv'),
  ('Gene-ranking agreement (rule-selected element)',"24",'rev24_gene_ranking_agreement.tsv'),
  ('GO-BP agreement (rule-selected element)',"24",'rev24_GO_agreement.tsv')]]))
S.append(('S7','Supplementary Table S7. E-P factors among the 100 largest core elements',
 ['Frequency of E-P factors, recurring (E-P, cell) combinations, cell types at the two ends of each cell factor, and gene rankings, GO-BP terms and enhancer-TSS distances of the most frequent E-P factors.',
  'Complete GO tables are written by Script 25 to the working directory and are not included here.'],
 [(n,src(25,f),rd(f)) for n,f in [
  ('E-P factor frequency','rev25_EP_factor_frequency.tsv'),
  ('Number of top-100 elements per (E-P, cell) combination','rev25_EP_by_CELL_counts.tsv'),
  ('Sum of |core| per (E-P, cell) combination','rev25_EP_by_CELL_abs_core.tsv'),
  ('Cell types at the two ends of each cell factor','rev25_cell_component_ends.tsv'),
  ('Top genes of the most frequent E-P factors','rev25_recurrent_EP_top_genes.tsv'),
  ('Top 10 GO-BP terms of the most frequent E-P factors','rev25_recurrent_EP_GO_top10.tsv'),
  ('Enhancer-TSS distance of the top-scoring pairs','rev25_recurrent_EP_distance.tsv')]]))
idx=wb.create_sheet('README',0)
idx['A1']='Supplementary Tables S1-S7'; idx['A1'].font=T
idx['A2']='Taguchi and Turki, Implicit hierarchical tensor decomposition of single-cell four-omics data (revised manuscript).'; idx['A2'].font=I
idx['A3']='Every table is an unmodified output of the R scripts in https://github.com/tagtag/TDCHARM (script number given with each table), produced by running them on GSE303006; each script also writes a log of its run.'; idx['A3'].font=I
idx.cell(5,1,'Sheet').font=B; idx.cell(5,2,'Content').font=B
for i,(n,t,d,b) in enumerate(S,6):
    idx.cell(i,1,n).font=F; idx.cell(i,2,t.split('. ',1)[1]).font=F
idx.column_dimensions['A'].width=10; idx.column_dimensions['B'].width=90
for n,t,d,b in S: sheet(n,t,d,b)
wb.save(OUT); print('saved',OUT,[ws.title for ws in wb])
