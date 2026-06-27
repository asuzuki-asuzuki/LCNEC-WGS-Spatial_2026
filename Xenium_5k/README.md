# Xenium 5k analysis

## Requirements
R (>=4.0; tested with 4.4.3)

The following R packages are required.
1. Seurat (v5.4.0)
2. Monocle 3 (v1.3.7)
3. CellChat (v2.2.0.9001)

It generally takes only several tens of minutes for installlation of these packages.  
For conducting analyses, it usually takes short time but sometimes it takes several hours depending on the number of samples. Memory requirements depend on the data size.

## Analysis
### Basic analysis of Xenium 5k data
- [Xenium_5k_Seurat.R](./Xenium_5k_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat in each sample.
- [Xenium_5k_Seurat_marker_basic.R](./Xenium_5k_Seurat_marker_basic.R): visualization of ASCL1/NEUROD1/POU2F3/YAP1 expression patterns by Seurat in each sample.
- [Xenium_5k_Seurat_merge.R](./Xenium_5k_Seurat_merge.R): merging Seurat objects of all samples.
- [Xenium_5k_Seurat_merge_sketch.R](./Xenium_5k_Seurat_merge_sketch.R): subsampling the cells and creating the sketch data.
- [Xenium_5k_Seurat_merge_sketch_ProjectData.R](./Xenium_5k_Seurat_merge_sketch_ProjectData.R): projection full merged data to the sketch data.
- [Xenium_5k_Seurat_merge_sketch_FindMarkers.R](./Xenium_5k_Seurat_merge_sketch_FindMarkers.R): extraction of DEGs in each cluster using the sketch data.
- [Xenium_5k_Seurat_merge_sketch_ProjectData_annotation.R](./Xenium_5k_Seurat_merge_sketch_ProjectData_annotation.R): set cell type annotation in each cluster.

### Trajectory analysis of epithelial cells
- [Xenium_5k_Seurat_sub.R](./Xenium_5k_Seurat_sub.R): extraction of epithelial sub-clusters by Seurat.
- [Monocle3_1_Xenium_5k.R](./Monocle3_1_Xenium_5k.R): consrtucting trajectories by Monocle 3.
- [Monocle3_2_Xenium_5k.R](./Monocle3_2_Xenium_5k.R): calculating pseudotime and generating plots.

### Cell-cell communication inference
- [CellChat_Xenium_5k](./CellChat_Xenium_5k.R): extraction of ligand-receptor interactions by CellChat.
- [CellChat_vis_Xenium_5k](./CellChat_vis_Xenium_5k.R): generating plots.

