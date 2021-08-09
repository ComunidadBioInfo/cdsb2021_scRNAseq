## ---- warning=FALSE, message=FALSE-----------------------------------------------------------
# Usemos datos de pbmc4k
library(BiocFileCache)
bfc <- BiocFileCache()
raw.path <- bfcrpath(bfc, file.path(
    "http://cf.10xgenomics.com/samples",
    "cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz"
))
untar(raw.path, exdir = file.path(tempdir(), "pbmc4k"))

library(DropletUtils)
library(Matrix)
fname <- file.path(tempdir(), "pbmc4k/raw_gene_bc_matrices/GRCh38")
sce.pbmc <- read10xCounts(fname, col.names = TRUE)
sce.pbmc


## ---- warning=FALSE, message=FALSE-----------------------------------------------------------
# Anotación de los genes
library(scater)
rownames(sce.pbmc) <- uniquifyFeatureNames(
    rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol
)
library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86,
    keys = rowData(sce.pbmc)$ID,
    column = "SEQNAME", keytype = "GENEID"
)

# Detección de _droplets_ con células
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]


## ---- warning=FALSE, message=FALSE-----------------------------------------------------------
# Control de calidad
stats <- perCellQCMetrics(sce.pbmc,
    subsets = list(Mito = which(location == "MT"))
)
high.mito <- isOutlier(stats$subsets_Mito_percent,
    type = "higher"
)
sce.pbmc <- sce.pbmc[, !high.mito]

# Normalización de los datos
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster = clusters)
sce.pbmc <- logNormCounts(sce.pbmc)


## ---- warning=FALSE, message=FALSE-----------------------------------------------------------
# Varianza de las log-counts
library(scran)
dec.pbmc <- modelGeneVar(sce.pbmc)


## ---- warning=FALSE, message=FALSE-----------------------------------------------------------
# Visualicemos la relación entre la media y la varianza
fit.pbmc <- metadata(dec.pbmc)
plot(fit.pbmc$mean, fit.pbmc$var,
    xlab = "Mean of log-expression",
    ylab = "Variance of log-expression"
)
curve(fit.pbmc$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)


## ---- warning=FALSE, message=FALSE-----------------------------------------------------------
# Ordenemos por los genes más interesantes para checar
# los datos
dec.pbmc[order(dec.pbmc$bio, decreasing = TRUE), ]


## ---- warning=FALSE, message=FALSE-----------------------------------------------------------
# Coeficiente de variación
dec.cv2.pbmc <- modelGeneCV2(sce.pbmc)


## ---- warning=FALSE, message=FALSE, echo = FALSE---------------------------------------------
# Visualicemos la relación con la media
fit.cv2.pbmc <- metadata(dec.cv2.pbmc)
plot(fit.cv2.pbmc$mean, fit.cv2.pbmc$cv2,
    log = "xy"
)
curve(fit.cv2.pbmc$trend(x),
    col = "dodgerblue",
    add = TRUE, lwd = 2
)


## ---- warning=FALSE, message=FALSE-----------------------------------------------------------
# Ordenemos por los genes más interesantes para checar
# los datos
dec.cv2.pbmc[order(dec.cv2.pbmc$ratio,
    decreasing = TRUE
), ]


## --------------------------------------------------------------------------------------------
## Información de la sesión de R
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

