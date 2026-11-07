# Implicit hierarchical tensor decomposition of GSE303006 single-cell four-omics

Y-h. Taguchi, Turki Turki

https://doi.org/10.64898/2026.08.18.745447


R code accompanying the manuscript:

> **Implicit Hierarchical Tensor Decomposition of Single-Cell Four-Omics Data Reveals Cell-Type-Linked Candidate Enhancer–Promoter-Associated Programs**

**Versions.** This branch holds the revised version (scripts 19–27 and `run_revision.R` added in response to the reviewers; git tag `v2-revision`). The code exactly as it accompanied the first bioRxiv version (v1) is preserved under the git tag `v1-preprint`.

The workflow analyzes the mouse-brain CHARM data set **GSE303006** using four
modalities measured in the same cells:

- ATAC
- H3K27ac
- RNA
- reconstructed 3D genome coordinates

The primary analysis uses a **non-K-mixed hierarchical HOSVD**. Each modality
is first reduced independently to 20 components. The resulting low-rank
representations are mapped to enhancer–promoter (E–P) pairs with
inverse-square-root degree correction and combined as an **implicit**
four-way tensor:

```text
E-P pair x cell x stage-1 component x modality
730,969 x 4,258 x 20 x 4
```

The full tensor is never materialized.

## Repository organization

```text
R/
  00_install_packages.R
  01_prepare_annotation_independent_inputs.R
  02_build_M23_EP_and_3D_summary.R
  03_finalize_EP_RNA_and_cell_QC.R
  04_build_3D_HDF5.R
  05_stage1_independent_K20.R
  06_sign_only_stage1_K20.R
  07_prepare_stage2_implicit.R
  08_stage2_implicit_HOSVD.R
  09_trace_core_to_original.R
  10_celltype_replicate_annotation.R
  11_core12_allEP_gene_rank.R
  12_core12_mouse_GO_BP_GSEA.R
  13_core12_GO_reduce_redundancy.R
  14_core12_GO_process_gene_overlap.R
  15_core12_original4modality_backprojection.R
  16_core12_replicate_reproducibility.R
  17_top100_CELLcomponent_effects.R
  18_top100_publication_figure.R
  revision_lib.R
  19_stage1_stage2_energy_report.R
  20_near_degenerate_rotation.R
  21_replicate_subsets.R
  21b_factor_matching.R
  21c_core12_by_position.R
  22_celltype_replicate_joint.R
  22b_contrast_by_replicate.R
  23_EP_link_benchmark.R
  23b_CHARM_link_concordance.R
  23c_CHARM_concordance_figure.R
  24_modality_ablation.R
  24b_ablation_factor_matching.R
  25_EP_factor_recurrence.R
  26_celltype_component_heatmap.R
  27_method_schematic.R

R_optional/
  core12_top50_target_gene_summary.R

run_preprocessing.R
run_analysis.R
run_revision.R

jobs/
  download_inputs.sh
  site_kyoto_icr.sh
  submit_all.sh
  tdcharm_step.sh

tests/
  make_synthetic_workdir.R

tools/
  build_supplementary_tables.py
```

The Procrustes-alignment experiment and the K-mixed HOSVD experiment explored
during method development are intentionally **not** part of the released
primary pipeline, because they were not used for the reported manuscript
results.

## Software requirements

Recommended:

- R >= 4.3
- Linux/Unix shell
- `zcat`
- `awk`
- MACS3 if ATAC/H3K27ac peaks are regenerated

Install/check R dependencies with:

```bash
Rscript R/00_install_packages.R
```

Main R dependencies include:

- data.table
- R.utils (used by data.table to read .gz files)
- Matrix
- RSpectra
- GenomicRanges
- rtracklayer
- rhdf5
- clusterProfiler
- org.Mm.eg.db
- AnnotationDbi
- GOSemSim
- GO.db
- igraph

## Expected input files

Download the corresponding GSE303006 supplementary files and place them in the
analysis working directory using these names:

```text
GSE303006_brain.atac.fragment.tsv.gz
GSE303006_brain.h3k27ac.fragment.tsv.gz
GSE303006_brain.RNAcounts.gene.total.format.tsv.gz
GSE303006_charm_metadata_qcpass.tsv.gz
gencode.vM23.annotation.gtf.gz
```

The reconstructed 20-kb 3D files are expected under:

```text
tdg20k/geo_submit/brain_20k_tdg/
```

with files named as `*.20k.3dg.gz`.

If peak calling is not rerun, the following peak files are also expected:

```text
brain_ATAC_peaks.narrowPeak
brain_H3K27ac_peaks.broadPeak
```

Script 01 contains the MACS3 commands used to regenerate these files. If the
peak files are absent (full recomputation from the GEO files), script 01 now
runs MACS3 itself. MACS3 can live in its own conda environment
(`conda create -n macs3 -c conda-forge -c bioconda python=3.12 macs3=3.1.0`); set
`TDCHARM_MACS3` to the full path of the executable, or let `SITE=kyoto_icr`
find `<conda base>/envs/macs3/bin/macs3` automatically. The MACS3 version, every command
that was run and the resulting peak counts are written to
`GSE303006_M23_MACS3_peak_calling_log.txt` in the working directory, and the
MACS3 console output to `brain_ATAC_macs3_callpeak.log` and
`brain_H3K27ac_macs3_callpeak.log`. Because peaks (and therefore the numbers
of enhancers and E-P pairs) can differ slightly between MACS3 versions, the
version recorded there should be reported in the Methods. Set
`TDCHARM_RUN_MACS3=1` to force peak calling even when the files exist.

`jobs/download_inputs.sh` downloads the GEO files listed above, the GENCODE
vM23 annotation and the reconstructed 20-kb 3D structures
(`GSE303006_brain.20k.tdg.tar.gz`, extracted to
`tdg20k/geo_submit/brain_20k_tdg/`) into `$GSE303006_WORKDIR`. Run it on a node
with internet access; files that are already complete are skipped, so it can
be re-run after an interruption.

## Final preprocessing thresholds

The final manuscript analysis used:

```text
RNA gene detected in >= 1% of original RNA cells
E-P genomic separation <= 1 Mb
3D bin-pair coverage >= 0.90
ATAC enhancer detection fraction >= 0.005
H3K27ac enhancer detection fraction >= 0.005

Cell QC:
ATAC enhancer-associated fragments >= 250
H3K27ac enhancer-associated fragments >= 100
RNA library size >= 500
RNA detected selected genes >= 500
```

The expected final dimensions are:

```text
730,969 E-P pairs
42,669 enhancers
16,239 genes
391,435 unique 20-kb bin pairs
4,258 cells
```

The generated 3D HDF5 file was approximately 5 GiB in the reported run.

## Running the pipeline

Set a separate working directory for downloaded data and generated large
intermediate files:

```bash
export GSE303006_WORKDIR=/path/to/GSE303006_work
```

### 1. Preprocessing

```bash
Rscript run_preprocessing.R
```

This creates the final E–P table, normalized RNA matrix, final cell-QC table,
and blockwise HDF5 3D matrix.

### 2. Hierarchical tensor analysis and downstream biology

```bash
Rscript run_analysis.R
```

This executes the manuscript analysis from Stage 1 through the publication
figure.

Because several steps are computationally expensive, individual scripts can
also be run separately:

```bash
cd "$GSE303006_WORKDIR"
Rscript /path/to/repository/R/08_stage2_implicit_HOSVD.R
```

## Key analysis choices

### Stage 1

Each modality is decomposed independently to `K = 20`.

- ATAC/H3K27ac: logCP10K -> feature-wise z score
- RNA: logCP10K -> feature-wise z score
- 3D: `-log(distance)` -> bin-pair-wise z score
- missing standardized 3D values are set to zero
- each retained rank-20 modality reconstruction is scaled to Frobenius norm 1

### Sign orientation

Only deterministic sign orientation is applied. No Procrustes rotation,
component mixing, or component permutation is used.

### Stage 2

The conceptual tensor is:

```text
730,969 x 4,258 x 20 x 4
```

For E–P pair `p`, cell `c`, component `k`, and modality `m`, degree correction
uses:

```text
ATAC/H3K27ac: 1 / sqrt(d_e)
RNA:          1 / sqrt(d_g)
3D:           1 / sqrt(d_b)
```

where `d_e`, `d_g`, and `d_b` are the numbers of final E–P rows sharing an
enhancer, gene, or 20-kb bin pair.

The primary HOSVD ranks are:

```text
20 x 20 x 20 x 4
```

The original `K=20` ordering is retained; the K mode is not mixed.

## Core-12 downstream analysis

The manuscript uses core rank 12 as a representative detailed example after
first showing that the top-100 core cell components are generally
cell-type-dominated rather than replicate-dominated.

The downstream workflow includes:

1. back-projection of large Tucker-core elements;
2. cell-type and replicate association;
3. degree-corrected all-gene ranking;
4. mouse GO Biological Process preranked GSEA;
5. GO semantic-similarity redundancy reduction;
6. shared leading-edge gene modules;
7. direct back-projection into original RNA/ATAC/H3K27ac/3D measurements;
8. replicate-wise reproducibility analysis;
9. top-100 CELL-component generalization figure.

## Running on a cluster with qsub

`jobs/submit_all.sh` submits the workflow as a chain of four jobs; each job
starts only after the previous one has finished successfully:

```text
synthetic      smoke test of 05-11, 17-18 and 19-27 on synthetic data
preprocessing  run_preprocessing.R  (01-04)
analysis       run_analysis.R       (05-18)
revision       run_revision.R       (19-27)
```

SGE / Univa-Altair Grid Engine, PBS Pro and Torque are supported; the
scheduler is detected automatically or set with `SCHEDULER=sge|pbspro|torque`.

```bash
# check the generated qsub commands first
DRY_RUN=1 GSE303006_WORKDIR=/path/to/workdir \
  TDCHARM_R_SETUP="module load R/4.3.3" bash jobs/submit_all.sh

# submit
GSE303006_WORKDIR=/path/to/workdir \
  TDCHARM_R_SETUP="module load R/4.3.3" bash jobs/submit_all.sh
```

Site profiles set suitable defaults for a given system. For the Kyoto
University ICR supercomputer (PBS Professional; SMALL queue for the short
synthetic test, APC queue for the real data, 8 cores / 128 GB / 336 h per job):

```bash
SITE=kyoto_icr DRY_RUN=1 GSE303006_WORKDIR=/path/to/workdir bash jobs/submit_all.sh
SITE=kyoto_icr GSE303006_WORKDIR=/path/to/workdir bash jobs/submit_all.sh
```

R is taken from the conda environment `r_env`: the profile looks up the conda
installation on the login node and runs
`source <conda>/etc/profile.d/conda.sh; conda activate r_env` inside each job
(another environment: `TDCHARM_CONDA_ENV=name`; conda not on PATH:
`CONDA_BASE=/path/to/miniconda3`). The R packages listed in
`R/00_install_packages.R` must be installed in that environment. Queue names and
limits can be checked with `qstat -Q`, `qstat -Qf APC` and `qstatmyjobs`.

Main settings (environment variables): `SITE`, `QUEUE` (or `QUEUE_<step>`), `NCPU` [4], `MEM` [64G],
`WALLTIME` [72:00:00], per-step overrides such as `MEM_revision=96G` or
`WALLTIME_analysis=120:00:00`, `STEPS` (e.g. `STEPS="revision"` to submit only
the revision step), `SYNTH_DIR` (synthetic test directory, outside the real
working directory), `SGE_PE` [smp] and `SGE_MEM_RESOURCE` [h_vmem] for Grid
Engine sites, and `EXTRA_QSUB_OPTS` for anything site-specific (project or
group codes). `TDCHARM_R_SETUP` is executed inside each job before R starts
(if `module` is not defined in batch shells, use e.g.
`"source /etc/profile.d/modules.sh; module load R/4.3.3"`). All other
`TDCHARM_*` variables set when calling the script are passed on to the jobs.

Job logs are written to `$GSE303006_WORKDIR/qsub_logs/`; each log starts with
the host, job ID, R version and the git state of the repository. On Grid
Engine a failed step exits with status 100 so that the remaining jobs of the
chain are held rather than started. `jobs/tdcharm_step.sh` can also be run
directly with bash (set `TDCHARM_STEP`, `TDCHARM_REPO`, `GSE303006_WORKDIR`).

A single script can be submitted on its own with `STEPS=script` and
`TDCHARM_SCRIPT=<file under R/>`, for example (after downloading the CHARM
Supplementary Tables workbook into the working directory):

```bash
SITE=kyoto_icr QUEUE=SMALL NCPU=4 MEM=32G WALLTIME=6:00:00 STEPS=script \
  TDCHARM_SCRIPT=23b_CHARM_link_concordance.R \
  GSE303006_WORKDIR=/path/to/workdir bash jobs/submit_all.sh
```

## Revision analyses

Scripts 19-27 were added in response to the reviewers' comments. They run
after the manuscript pipeline, in the same working directory:

```bash
Rscript run_preprocessing.R
Rscript run_analysis.R
Rscript run_revision.R
```

Each script writes its tables and figures with the prefix
`GSE303006_M23_rev<NN>_` and saves its full console output as
`GSE303006_M23_rev<NN>_..._log.txt` (or `..._log.txt`) next to them, including
`sessionInfo()`. Shared functions are in `R/revision_lib.R`; its Stage-2
implementation follows the same algebra as script 08, and scripts 20 and 24
check it against the saved output of 08 (the check is printed in their logs).

| Script | Reviewer comment | What it computes |
|---|---|---|
| `19_stage1_stage2_energy_report.R` | R1-6, R3-1 | Stage-1 variance retention and reconstruction error per modality (with cumulative explained fraction for k = 1..20); Stage-2 retained fraction and reconstruction error; Stage-2 and overall retention per modality; check that the degree-corrected mapping preserves each modality's norm. No decomposition is recomputed. |
| `20_near_degenerate_rotation.R` | R1-2 | Blocks of nearly degenerate Stage-1 singular values (relative gap < 0.05); random (and 45-degree) rotations inside the blocks; Stage-2 recomputed; subspace angles, matched core and agreement of the core-#12 gene ranking. |
| `21_replicate_subsets.R` | R1-3 | Scripts 05-08 re-run unchanged on each two-replicate subset (R1+R2, R1+R3, R2+R3) in sub-directories `rev21_*`; the same pre-specified core-selection rule in every run; gene rankings, GO-BP GSEA, cell-side statistics, and their agreement among the subsets and with the full data. |
| `21b_factor_matching.R` | R1-3 | Direct comparison of the Stage-2 factor vectors of the full data and of each two-replicate subset (|correlation| of every E-P factor over all E-P pairs and of every cell factor over shared cells), independent of factor indices and of any selection rule. |
| `21c_core12_by_position.R` | R1-3 | The core element at the position of manuscript core #12 (EP03, CELL03, K02, MOD01; each run's own factor ranks) analysed in the full data and in every subset: core value and its rank, gene rankings, GO-BP GSEA, and their agreement. |
| `22_celltype_replicate_joint.R` | R1-4 | Cell-type x replicate composition; joint linear model (cell type, replicate, interaction, depth/QC covariates) with type-II partial eta^2 for every cell component (raw and rank-transformed loadings); within-cell-type replicate Kruskal-Wallis tests. |
| `22b_contrast_by_replicate.R` | R1-4 | For the cell component of core #12 and its two contrasted cell types: loadings per cell type and replicate, the contrast within each replicate, pairwise replicate comparisons within each cell type, and the largest replicate shift relative to the between-type difference. |
| `23_EP_link_benchmark.R` | R1-5 | Enrichment of core-#12-prioritized E-P pairs among supported links (ABC-like links from this data set; optional external link tables via `TDCHARM_EXTERNAL_LINKS`), stratified by enhancer-TSS distance (CMH test, AUPRC against a distance-matched permutation null). |
| `23b_CHARM_link_concordance.R` | R1-5 | Concordance of core-#12-prioritized E-P pairs with the E-P links of the original CHARM study (Chen et al., Nature 2026, Supplementary Table 7; links in any cell type, intergenic only, and in the two contrasted cell types), in all candidate pairs and in pairs of genes linked by CHARM; same CMH/AUPRC procedures as 23, plus core-#12 scores of Ndnf/Lamp5-specific versus L3/4 IT-specific CHARM links. The workbook `41586_2026_10322_MOESM3_ESM.xlsx` is not redistributed: download it from the article page into the working directory (or set `TDCHARM_CHARM_SUPP`); otherwise the step is skipped. These links are predictions from the same data set, so this is a concordance analysis, not a validation. |
| `23c_CHARM_concordance_figure.R` | R1-5 | Figure from the outputs of 23b: (a) core-#12 score percentile of candidate pairs linked by CHARM only in Ndnf/Lamp5, only in L3/4 IT, or only in other cell types; (b) Ndnf/Lamp5 and L3/4 IT links among the top 1,000 pairs on each side of the ranking, with expected counts and distance-stratified CMH odds ratios. |
| `24_modality_ablation.R` | R1-1 | Stage 2 recomputed for each single modality and each leave-one-modality-out set; the same core-selection rule; agreement of gene rankings and GO with the four-modality result. |
| `24b_ablation_factor_matching.R` | R1-1 | Direct comparison of the factor vectors of every modality variant of 24 with the four-modality factors (|correlation| of every E-P factor over all E-P pairs, accumulated block-wise, and of every cell factor), independent of any selection rule; the core element built from the factors matching core #12 in each variant (value, rank, gene ranking, GO-BP GSEA) and its agreement with core #12. Modality columns are sign-aligned to the core-#12 modality column over the shared modalities. |
| `25_EP_factor_recurrence.R` | R2-4 | Frequency of E-P factors among the top-100 core elements; recurring (E-P, cell) factor combinations; gene ranking, GO and E-P distance for the most frequent E-P factors. |
| `26_celltype_component_heatmap.R` | R2-2 | Cell type x cell component heatmap (mean standardized loading) for all cell components in the top-100 core elements. |
| `27_method_schematic.R` | R2 minor 1 | Graphical summary of the method (draft schematic). |

Figure revisions in the original scripts:

- `15_...backprojection.R`: module-score heatmap (manuscript Fig. 3) uses a
  higher-contrast diverging palette with a colour bar and labelled cell groups
  (R3-3); also saved as PNG.
- `16_...replicate_reproducibility.R`: replicate effect heatmap (Fig. 4) has a
  colour bar and the effect values printed in the cells (R2 minor).
- `18_top100_publication_figure.R`: the region label that overlapped the
  sign-test summary (the "shadow") was moved (R2 minor).

Pre-specified core-selection rule (scripts 21 and 24). In each run, the cell
component whose 50 most positive and 50 most negative cells are most enriched
for the two cell types that characterize core #12 (Inh Ndnf/Lamp5 and Ex L3/4
IT) is selected, then the largest |core| element containing it. The labels are
matched by regular expression against the cell-type column chosen by script 10;
set `TDCHARM_GROUP_A` / `TDCHARM_GROUP_B` to exact labels if needed. Script 21
reports whether the rule reproduces manuscript core #12 on the full data.

Optional settings: `TDCHARM_RUN_GSEA=0` skips GSEA (21, 24, 25);
`TDCHARM_RECOMPUTE=1` recomputes the cached subset runs of 21;
`TDCHARM_N_ROT`, `TDCHARM_GAP_TOL` (20); `TDCHARM_N_PERM`, `TDCHARM_TOP_N`,
`TDCHARM_ABC_THRESHOLD` (23); `TDCHARM_N_RECURRENT_EP` (25).

### Supplementary Tables S1-S7

After `run_revision.R`, the tables of the revised manuscript are assembled
into one workbook (one sheet per table, each labelled with its source script
and file; requires Python 3 with pandas and openpyxl):

```bash
python3 tools/build_supplementary_tables.py /path/to/workdir Supplementary_Tables_S1-S7.xlsx
```

Large outputs (per-gene scores, complete GO tables, factor correlation
matrices, ABC-like scores of all pairs) stay in the working directory.

### Synthetic test data

`tests/make_synthetic_workdir.R` writes a small synthetic working directory
with the same file names and object structures as the real intermediate
files. It was used to check that scripts 05-11, 17-18 and 19-27 run end to
end, that the shared Stage-2 implementation reproduces script 08 exactly,
and that the shared gene score reproduces script 11. Numbers obtained from
the synthetic data have no biological meaning.

```bash
Rscript tests/make_synthetic_workdir.R /tmp/tdcharm_synthetic
```

## Reproducibility notes

The R scripts operate on large intermediate files and therefore do not ship
data in the Git repository. `.gitignore` excludes RDS/HDF5/TSV outputs,
downloaded raw data, and generated figures.

The R1/R2/R3 analysis is a **within-dataset reproducibility analysis**, not
an independent validation, because the HOSVD factors were learned using the
combined data set.

## Code history

The released scripts correspond to the final manuscript analysis. Development
experiments that were ultimately rejected from the primary method
(Procrustes component alignment and K-mixed stage-2 HOSVD) are not included
in the main workflow.

## Citation

If you use this workflow, please cite the associated manuscript once the
final bibliographic information is available.

## License

No software license is assigned in this draft repository. Add the license
chosen by the author before public release.
