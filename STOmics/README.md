# Visium analysis

## Requirements
R (>=4.0; tested with 4.3.2)  

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
