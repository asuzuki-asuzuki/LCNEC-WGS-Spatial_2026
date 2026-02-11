# STOmics analysis

## Requirements
R (>=4.0; tested with 4.3.2)  

The following R packages are required.
1. Seurat (v5.0.3)
  
It generally takes only several tens of minutes for installlation of these packages.  
For conducting analyses, it usually takes short time (up to several hours). Memory requirements depend on the data size.  

## Analysis
### Basic analysis of STOmics data (bin 50)
- [STOmics_bin50_Seurat.R](./STOmics_bin50_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat.
- [STOmics_bin50_Seurat_marker.R](./STOmics_bin50_Seurat_marker.R): generating spatial plots with marker gene expression.
- [STOmics_bin50_Seurat_Module_projection.R](./STOmics_bin50_Seurat_Module_projection.R): generating UMAP and spatial plots with co-expression network module scores.
