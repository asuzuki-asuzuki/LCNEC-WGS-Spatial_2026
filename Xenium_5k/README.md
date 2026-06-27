# Xenium 5k analysis

## Requirements
R (>=4.0; tested with 4.4.3)

The following R packages are required.
1. Seurat (v5.4.0)
2. Monocle 3 (v1.3.7)

It generally takes only several tens of minutes for installlation of these packages.  
For conducting analyses, it usually takes short time but sometimes it takes several hours depending on the number of samples. Memory requirements depend on the data size.

## Analysis
### Basic analysis of Xenium data
- [Xenium_k_Seurat.R](./Xenium_5k_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat in each sample.
- [Xenium_5k_Seurat_marker_basic.R](./Xenium_5k_Seurat_marker_basic.R): visualization of ASCL1/NEUROD1/POU2F3/YAP1 expression patterns by Seurat in each sample.
- [Xenium_5k_Seurat_merge.R](./Xenium_5k_Seurat_merge.R): merging Seurat objects of all samples.
