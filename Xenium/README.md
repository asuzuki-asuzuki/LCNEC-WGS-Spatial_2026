# Xenium v1 analysis

## Requirements
R (>=4.0; tested with 4.3.2)

The following R packages are required.
1. Seurat (v5.0.3)
2. Monocle 3 (v1.3.7)

It generally takes only several tens of minutes for installlation of these packages.  
For conducting analyses, it usually takes short time but sometimes it takes several hours depending on the number of samples. Memory requirements depend on the data size.

## Analysis
### Basic analysis of Xenium data
- [Xenium_Seurat.R](./Xenium_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat in each sample.
- [Xenium_Seurat_merge.R](./Xenium_Seurat_merge.R): merging Seurat objects of all samples.
- [Xenium_Seurat_merge_sketch.R](./Xenium_Seurat_merge_sketch.R): basic analysis for subsampled data (sketch analysis).
- [Xenium_Seurat_merge_sketch_ProjectData.R](./Xenium_Seurat_merge_sketch_ProjectData.R): projection full merged data to the sketch data.
- [Xenium_Seurat_merge_sketch_FindMarkers.R](./Xenium_Seurat_merge_sketch_FindMarkers.R): extraction of DEGs in each cluster using the sketch data.
- [Xenium_Seurat_merge_featureplot.R](./Xenium_Seurat_merge_featureplot.R): maker gene expression (raw count) on the UMAP plot in the sketch and full data.
- [Xenium_Seurat_merge_sketch_ProjectData_annotation.R](./Xenium_Seurat_merge_sketch_ProjectData_annotation.R): set cell type annotation in each cluster.

### Analysis of epithelial cells
- [Xenium_Seurat_sub_sketch.R](./Xenium_Seurat_sub_sketch.R): basic analysis of epithelial sub-clusters by Seurat.
- [Xenium_Seurat_sub_sketch_FindMarkers.R](./Xenium_Seurat_sub_sketch_FindMarkers.R): extraction of DEGs in each cluster using the sketch data.
- [Xenium_Seurat_sub_sketch_ProjectData_subannotation.R](./Xenium_Seurat_sub_sketch_ProjectData_subannotation.R): projection full merged data to the sketch data and set cell subtype annotation in each cluster.
- [Xenium_Seurat_VlnPlot_limit.R](./Xenium_Seurat_VlnPlot_limit.R): ASCL1, ASCL2 and YAP1 expression on violin plots in each cell sub-clusters.
- [Xenium_Seurat_count_ASCL1_ASCL2_YAP1.R](./Xenium_Seurat_count_ASCL1_ASCL2_YAP1.R): counting ASCL1, ASCL2 and YAP1-positive cells in tumor cells.

### Trajectory analysis of epithelial cells
- [Monocle3_1.R](./Monocle3_1.R): consrtucting trajectories by Monocle 3.
- [Monocle3_2.R](./Monocle3_2.R): calculating pseudotime and generating plots.
- [Monocle3_3.R](./Monocle3_3.R): generating plots.

### Cellular niche analysis
<I>Preproceccing and Niche analysis for each sample</I><br>
- [Xenium_Seurat_extract_ID2annotation_subannotation.R](./Xenium_Seurat_extract_ID2annotation_subannotation.R): extraction of cell type annotation and subannotation in each cell.
- [Xenium_Seurat_extract_ID2annotation_subannotation_list.pl](./Xenium_Seurat_extract_ID2annotation_subannotation.pl): extraction of cell type annotation and subannotation in each cell for each sample.
- [Xenium_Seurat_annotation.R](./Xenium_Seurat_annotation.R): setting annotation and subannotation of each cell in each sample's object.
- [niche_TME.R](./niche_TME.R): niche analysis for each sample.

<I>Niche analysis for all samples</I><br>
- [Xenium_Seurat_merge_TME.R](./Xenium_Seurat_merge_TME.R): merging Seurat objects with the niche matrix of all samples.
- [niche_merge_TME.R](./niche_nerge_TME.R): niche analysis for all samples.









