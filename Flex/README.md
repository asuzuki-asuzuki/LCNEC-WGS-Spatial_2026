# ASCL1 siRNA knockdown scRNA-seq analysis

## Requirements
R (>=4.0; tested with 4.4.3)  

The following R packages are required.
1. Seurat (v5.3.1)
2. Slingshot (v2.14.0)
3. tradeSeq (v1.20.0) 

It generally takes only several tens of minutes for installation of these packages.  
For conducting analyses, it usually takes short time (up to several hours). Memory requirements depend on the data size.  

## Analysis
### Basic analysis of scRNA-seq data
- [01_qc.R](./01_qc.R): checking QC metrics (nCount, nFeature, and percent.mt).
- [02_filter_umap.R](./02_filter_umap.R): filtering, normalization, regression of cell-cycle scores, dimensional reduction, clustering, and DEG extraction.
- [03_marker_plot.R](./03_marker_plot.R): visualization of marker genes and module scores.
- [04_subset.R](./04_subset.R): subsetting of selected clusters and recomputing UMAP.

### Trajectory and pseudotime analysis of cell lines
- [05_trajectory.R](./05_trajectory.R): constructing trajectories and calculating pseudotime by Slingshot.
- [06_tradeseq_pseudotime.R](./06_tradeseq_pseudotime.R): identification of pseudotime-associated genes by tradeSeq.
- [07_pseudotime_modules.R](./07_pseudotime_modules.R): clustering of pseudotime-associated genes into gene modules.

### Mapping gene modules onto tissue data
- [08_module_score_H810.R](./08_module_score_H810.R): scoring of H810-derived gene modules on the Xenium 5k UMAP.
- [08_module_score_VMRC.R](./08_module_score_VMRC.R): scoring of VMRC-LCD-derived gene modules on the Xenium 5k UMAP.
