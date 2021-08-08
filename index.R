## ----install, eval = FALSE-------------------------------------------------------------------------------------------
## ## Para instalar paquetes
## if (!requireNamespace("remotes", quietly = TRUE)) {
##     install.packages("remotes")
## }
## 
## ## Para instalar paquetes de Bioconductor
## remotes::install_cran("BiocManager")
## BiocManager::version()
## # El anterior comando debe mostrar que estás usando la versión 3.13
## 
## ## Instala los paquetes de R que necesitamos
## BiocManager::install(
##     c(
##         "SingleCellExperiment",
##         "usethis",
##         "here",
##         "scran",
##         "scater",
##         "scRNAseq",
##         "org.Mm.eg.db",
##         "AnnotationHub",
##         "ExperimentHub",
##         "BiocFileCache",
##         "DropletUtils",
##         "EnsDb.Hsapiens.v86",
##         "TENxPBMCData",
##         "BiocSingular",
##         "batchelor",
##         "uwot",
##         "Rtsne",
##         "pheatmap",
##         "fossil",
##         "ggplot2",
##         "cowplot",
##         "RColorBrewer",
##         "plotly",
##         "iSEE",
##         "pryr",
##         "sessioninfo",
##         "scPipe",
##         "Seurat"
##     )
## )


## ----session_packages, eval = TRUE, message = FALSE------------------------------------------------------------------
library("SingleCellExperiment")
library("usethis")
library("here")
library("scran")
library("scater")
library("scRNAseq")
library("org.Mm.eg.db")
library("AnnotationHub")
library("ExperimentHub")
library("BiocFileCache")
library("DropletUtils")
library("EnsDb.Hsapiens.v86")
library("TENxPBMCData")
library("BiocSingular")
library("batchelor")
library("uwot")
library("Rtsne")
library("pheatmap")
library("fossil")
library("ggplot2")
library("cowplot")
library("RColorBrewer")
library("plotly")
library("iSEE")
library("pryr")
library("sessioninfo")
library("scPipe")
library("PCAtools")
library("bluster")


## ----session_info----------------------------------------------------------------------------------------------------
## Reproducibility information
options(width = 120)
session_info()
proc.time()

