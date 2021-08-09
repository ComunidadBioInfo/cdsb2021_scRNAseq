## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
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


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
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


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
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


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
## Identificación de genes altamente variables
set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop = 0.1)


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
## Reducción de dimensiones
set.seed(10000)
sce.pbmc <- denoisePCA(sce.pbmc,
    subset.row = top.pbmc,
    technical = dec.pbmc
)
set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred = "PCA")
set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred = "PCA")


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
# clustering
g <- buildSNNGraph(sce.pbmc, k = 10, use.dimred = "PCA")
clust <- igraph::cluster_walktrap(g)$membership
sce.pbmc$cluster <- factor(clust)


## ---- warning=FALSE, message=FALSE, eval = FALSE-----------------------------------------------------
## # Human
## SingleR::BlueprintEncodeData()
## SingleR::DatabaseImmuneCellExpressionData()
## SingleR::HumanPrimaryCellAtlasData()
## SingleR::MonacoImmuneData()
## SingleR::NovershternHematopoieticData()
## 
## # Mice
## SingleR::ImmGenData()
## SingleR::MouseRNASeqData()


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
# if needed install celldex
# create directory? y
library(celldex)
ref <- celldex::BlueprintEncodeData()


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
library(SingleR)
pred <- SingleR(
    test = sce.pbmc, ref = ref,
    labels = ref$label.main
)


## ---- warning=FALSE, message=FALSE, echo=FALSE-------------------------------------------------------
plotScoreHeatmap(pred)


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
total_pruned <- sum(is.na(pred$pruned.labels))
plotScoreHeatmap(pred, show.pruned = TRUE)


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
plotScoreDistribution(pred)


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
# install gmp, ClusterR, mbkmeans dependencies if needed
library(clusterExperiment)
sce.pbmc$labels <- pred$labels
all.markers <- metadata(pred)$de.genes
lab <- "B-cells"
# Get top-10 marker genes for B-cells compared to each other cell
# type
top.markers <- Reduce(union, sapply(all.markers[[lab]], head, 10))

# plotHeatmap(sce.pbmc, order_columns_by="labels",
#  features=top.markers, center=TRUE, zlim=c(-3, 3), main=lab)


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
tab <- table(Assigned = pred$pruned.labels, Cluster = sce.pbmc$cluster)

library(pheatmap)
# Proportion of cells in each cluster assigned to each label
pheatmap(prop.table(tab, margin = 2),
    color = colorRampPalette(c("white", "blue"))(101)
)
# (log-)number of cells in each cluster assigned to each label
# Adding a pseudo-count of 10 to avoid strong color jumps with just
# 1 cell.
pheatmap(log2(tab + 10),
    color = colorRampPalette(c("white", "blue"))(101)
)


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
plotTSNE(sce.pbmc, colour_by = "labels", text_by = "labels")


## ----------------------------------------------------------------------------------------------------
## Información de la sesión de R
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

