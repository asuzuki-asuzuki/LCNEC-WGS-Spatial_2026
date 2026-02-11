# Visium analysis

## Requirements
R (>=4.0; tested with 4.3.2)  
Python (3.7 or later)  

The following R packages are required.
1. Seurat (v5.0.3)
2. hdf5r (for loding HDF5 files)
3. hdWGCNA (v0.3.3)
4. Slingshot (v2.10.0)
  
It generally takes only several tens of minutes for installlation of these packages.  
For conducting analyses, it usually takes short time (up to several hours). Memory requirements depend on the data size.  

## Analysis
### Basic analysis of Visium data  
- [Visium_Seurat.R](./Visium_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat.
- [Visium_Seurat_merge.R](./Visium_Seurat_merge.R): merging Seurat objects of all samples.
- [Visium_Seurat_merge_harmony.R](./Visium_Seurat_merge_harmony.R): Harmony integration of all samples.
- [Visium_Seurat_merge_harmony_FindMarkers.R](./Visium_Seurat_merge_harmony_FindMarkers.R): extraction of DEGs in each cluster using the harmony integrated data.
- [Visium_Seurat_merge_harmony_annotation.R](./Visium_Seurat_merge_harmony_annotation.R): set cell type annotation in each cluster.

### Co-expression network analysis of epithelial cells
- [Visium_Seurat_sub.R](./Visium_Seurat_sub.R): basic analysis of epithelial sub-clusters by Seurat.
- [hdWGCNA.R](./hdWGCNA.R): construction of co-expression networks using epithelial sub-clusters by hdWGCNA.
- [hdWGCNA_network_visualization.R](./hdWGCNA_network_visualization.R): visualization of co-expression networks.
- [hdWGCNA_spatialplot.R](./hdWGCNA_spatialplot.R): visualization of module eigengenes in spatial plots.
- [hdWGCNA_enrichment.R](./hdWGCNA_enrichment.R): GO enrichment analysis of module members.
- [hdWGCNA_DME_DP_vs_Others.R](./hdWGCNA_DME_DP_vs_Others.R): comparison of module eigengenes between DP and other subtypes.

### Trajectory analysis of tumor cells
- [Visium_Seurat_sub_tumor.R](./Visium_Seurat_sub_tumor.R): basic analysis of tumor sub-clusters by Seurat.
- [Slingshot_1.R](./Slingshot_1.R): consrtucting trajectories by Slingshot.
- [Slingshot_2.R](./Slingshot_2.R): calculating pseudotime and generating plots.
- [Slingshot_3.R](./Slingshot_3.R): generating plots for a sample.



