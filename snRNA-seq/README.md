# snRNA-seq analysis

## Requirements
R (>=4.0; tested with 4.3.2)

The following R packages are required.
1. Seurat (v5.0.3)
2. Monocle 3 (v1.3.7)

It generally takes only several tens of minutes for installlation of these packages.  
For conducting analyses, it usually takes short time but sometimes it takes several hours depending on the number of samples. Memory requirements depend on the data size.

## Analysis
### Basic analysis of snRNA-seq data in each sample
- [scRNA_Seurat_qc.R](./scRNA_Seurat_qc.R): generation of QC plots by Seurat in each sample.
- [scRNA_Seurat.R](./scRNA_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat in each sample.
- [scRNA_Seurat_annotation.R](./scRNA_Seurat_annotation.R): set cell type annotation in each cluster.

### Basic analysis of snRNA-seq data in all samples
- [scRNA_Seurat_merge.R](./scRNA_Seurat_merge.R): merging Seurat objects of all samples.
- [scRNA_Seurat_merge_harmony.R](./scRNA_Seurat_harmony.R): Harmony integration of all samples.
- [scRNA_Seurat_marker_ADAMTS.R](./scRNA_Seurat_marker_ADAMTS.R): Visualization of ADAMTS family gene expressions.

### Trajectory analysis
- [Monocle3_1_all.R](./Monocle3_1_all.R): consrtucting trajectories by Monocle 3.
- [Monocle3_2_all.R](./Monocle3_2_all.R): calculating pseudotime and generating plots.
