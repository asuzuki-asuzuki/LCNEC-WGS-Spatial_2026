# Visium analysis

## Requirements
R (>=4.0; tested with 4.0.2 and 4.2.1)  
Python (3.7 or later)  

The following R packages are required.
1. Seurat (v4)
2. hdf5r (for loding HDF5 files)
  
It generally takes only several tens of minutes for installlation of these packages.  
For conducting analyses, it usually takes short time (up to several hours). Memory requirements depend on the data size.  

## Analysis
### Basic analysis of Visium data  
- [Visium_Seurat.R](./Visium_Seurat.R): basic analysis (dimentional reduction, clustering, etc.) by Seurat.

### Visualization

