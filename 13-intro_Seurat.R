## ---- warning=FALSE, message=FALSE------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
# Load the PBMC dataset
proydir <- "/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/CDSB/clustering/"
pbmc.data <- Read10X(data.dir = paste0(proydir, "data/filtered_gene_bc_matrices/hg19/"))
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc


## ---- warning=FALSE, message=FALSE------------------------------------------------------------
str(pbmc)


## ---- warning=FALSE, message=FALSE------------------------------------------------------------
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size
dense.size / sparse.size


## ---------------------------------------------------------------------------------------------
## Información de la sesión de R
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

