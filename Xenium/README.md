# Xenium v1 analysis

## Requirements
R (>=4.0; tested with 4.3.2)

The following R packages are required.
1. Seurat (v5.0.3)

It generally takes only several tens of minutes for installlation of these packages.  
For conducting analyses, it usually takes short time but sometimes it takes several hours depending on the number of samples. Memory requirements depend on the data size.

## Analysis
### Basic analysis of Xenium data
- [Xenium_Seurat.R](./Xenium_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat in each sample.
- [Xenium_Seurat_merge.R](./Xenium_Seurat_merge.R): merging Seurat objects of all samples
- [Xenium_Seurat_merge_sketch.R](./Xenium_Seurat_merge_sketch.R): basic analysis for subsampled data (sketch analysis) 
- [Xenium_Seurat_merge_sketch_ProjectData.R](./Xenium_Seurat_merge_sketch_ProjectData.R): projection full merged data to the sketch data
- [Xenium_Seurat_merge_sketch_FindMarkers.R](./Xenium_Seurat_merge_sketch_FindMarkers.R): extraction of DEGs in each cluster using the sketch data
- [Xenium_Seurat_merge_featureplot.R](./Xenium_Seurat_merge_featureplot.R): maker gene expression (raw count) on the UMAP plot in the sketch and full data
- [Xenium_Seurat_merge_sketch_ProjectData_annotation.R](./Xenium_Seurat_merge_sketch_ProjectData_annotation.R): set cell type annotation in each cluster

### Analysis of epithelial cells
- [Xenium_Seurat.R](./Xenium_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat in each sample.
