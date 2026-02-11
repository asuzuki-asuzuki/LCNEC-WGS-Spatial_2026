# Visium analysis

## Requirements
R (>=4.0; tested with 4.3.2)  
Python (3.7 or later)  

The following R packages are required.
1. Seurat (v5.0.3)
2. hdf5r (for loding HDF5 files)
  
It generally takes only several tens of minutes for installlation of these packages.  
For conducting analyses, it usually takes short time (up to several hours). Memory requirements depend on the data size.  

## Analysis
### Basic analysis of Visium data  
- [Visium_Seurat.R](./Visium_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat.
- [Visium_Seurat_merge.R](./Visium_Seurat_merge.R): merging Seurat objects of all samples.
- [Visium_Seurat_merge_harmony.R](./Visium_Seurat_merge_harmony.R): Harmony integration of all samples.
- [Visium_Seurat_merge_harmony_FindMarkers.R](./Visium_Seurat_merge_harmony_FindMarkers.R): extraction of DEGs in each cluster using the harmony integrated data.


### Visualization

