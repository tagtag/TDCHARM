# Implicit hierarchical tensor decomposition of GSE303006 single-cell four-omics

R code accompanying the manuscript:

> **Implicit Hierarchical Tensor Decomposition of Single-Cell Four-Omics Data Reveals Cell-Type-Associated Enhancer–Promoter Regulatory Programs**

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

R_optional/
  core12_top50_target_gene_summary.R

run_preprocessing.R
run_analysis.R
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

Script 01 contains the MACS3 commands used to regenerate these files.

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
