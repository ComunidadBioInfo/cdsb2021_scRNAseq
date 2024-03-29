## ---- warning=FALSE, message=FALSE----------------------
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


## ---- warning=FALSE, message=FALSE----------------------
# gene-annotation
library(scater)
rownames(sce.pbmc) <- uniquifyFeatureNames(
    rowData(sce.pbmc)$ID, rowData(sce.pbmc)$Symbol
)
library(EnsDb.Hsapiens.v86)
location <- mapIds(EnsDb.Hsapiens.v86,
    keys = rowData(sce.pbmc)$ID,
    column = "SEQNAME", keytype = "GENEID"
)

# cell-detection
set.seed(100)
e.out <- emptyDrops(counts(sce.pbmc))
sce.pbmc <- sce.pbmc[, which(e.out$FDR <= 0.001)]


## ---- warning=FALSE, message=FALSE----------------------
# quality-control
stats <- perCellQCMetrics(sce.pbmc,
    subsets = list(Mito = which(location == "MT"))
)
high.mito <- isOutlier(stats$subsets_Mito_percent,
    type = "higher"
)
sce.pbmc <- sce.pbmc[, !high.mito]

# normalization
library(scran)
set.seed(1000)
clusters <- quickCluster(sce.pbmc)
sce.pbmc <- computeSumFactors(sce.pbmc, cluster = clusters)
sce.pbmc <- logNormCounts(sce.pbmc)


## ---- warning=FALSE, message=FALSE----------------------
# variance modelling
set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop = 0.1)


## ---- warning=FALSE, message=FALSE----------------------
# dimensionality-reduction
set.seed(10000)
sce.pbmc <- denoisePCA(sce.pbmc,
    subset.row = top.pbmc,
    technical = dec.pbmc
)

set.seed(100000)
sce.pbmc <- runTSNE(sce.pbmc, dimred = "PCA")

set.seed(1000000)
sce.pbmc <- runUMAP(sce.pbmc, dimred = "PCA")


## ---- warning=FALSE, message=FALSE----------------------
library(scran)
# Build graph using k = 10 nearest neighbours in PCA-space
g <- buildSNNGraph(sce.pbmc, k = 10, use.dimred = "PCA")
# Identify communities using the Walktrap method
clust <- igraph::cluster_walktrap(g)$membership


## ---- warning=FALSE, message=FALSE, , fig.dim = c(6, 4)----
# Visualise clusters on t-SNE plot
library(scater)
sce.pbmc$cluster <- factor(clust)
plotReducedDim(sce.pbmc, "TSNE", colour_by = "cluster")


## ---- warning=FALSE, message=FALSE----------------------
library(scran)
# Build graph using k = 50 nearest neighbours in PCA-space
g50 <- buildSNNGraph(sce.pbmc, k = 50, use.dimred = "PCA")
# Identify communities using the Walktrap method
clust50 <- igraph::cluster_walktrap(g50)$membership


## ---- warning=FALSE, message=FALSE, , fig.dim = c(6, 4)----
# Visualise clusters on t-SNE plot
library(scater)
sce.pbmc$cluster50 <- factor(clust50)
plotReducedDim(sce.pbmc, "TSNE", colour_by = "cluster50")


## ---- warning=FALSE, message=FALSE , fig.dim = c(6, 4)----
# Jaccard-based weights followed by Louvain clustering
# aka 'Seurat-style' clustering
g2 <- buildSNNGraph(sce.pbmc, k = 10, use.dimred = "PCA", type = "jaccard")
clust2 <- igraph::cluster_louvain(g2)$membership
sce.pbmc$cluster2 <- factor(clust2)
plotReducedDim(sce.pbmc, "TSNE", colour_by = "cluster2")


## ---- fig.cap = "Estilo scran vs estilo Seurat.", fig.width=10----
library("patchwork")

plotReducedDim(sce.pbmc, "TSNE", colour_by = "cluster") +
    plotReducedDim(sce.pbmc, "TSNE", colour_by = "cluster2")


## ---- eval = FALSE--------------------------------------
## g.num <- buildSNNGraph(sce.pbmc, use.dimred = "PCA", type = "number")
## g.jaccard <- buildSNNGraph(sce.pbmc, use.dimred = "PCA", type = "jaccard")
## g.none <- buildKNNGraph(sce.pbmc, use.dimred = "PCA")


## ---- eval = FALSE--------------------------------------
## clust.louvain <- igraph::cluster_louvain(g)$membership
## clust.infomap <- igraph::cluster_infomap(g)$membership
## clust.fast <- igraph::cluster_fast_greedy(g)$membership
## clust.labprop <- igraph::cluster_label_prop(g)$membership
## clust.eigen <- igraph::cluster_leading_eigen(g)$membership


## ---- warning=FALSE, message=FALSE----------------------
library(bluster)

# obteniendo la métrica de modularidad
ratio <- pairwiseModularity(g, clust, as.ratio = TRUE)
dim(ratio)


## ---- warning=FALSE, message=FALSE, fig.dim = c(6, 4)----
library(pheatmap)
pheatmap(log2(ratio + 1),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", "blue"))(100)
)


## ---- warning=FALSE, message=FALSE----------------------
myClusterFUN <- function(x) {
    g <- buildSNNGraph(x, use.dimred = "PCA", type = "jaccard")
    igraph::cluster_louvain(g)$membership
}

originals <- myClusterFUN(sce.pbmc)
set.seed(0010010100)
coassign <- bootstrapStability(sce.pbmc,
    FUN = myClusterFUN,
    clusters = originals
)


## ---- warning=FALSE, message=FALSE , fig.dim = c(6, 4)----
pheatmap(coassign,
    cluster_row = FALSE, cluster_col = FALSE,
    color = rev(viridis::magma(100))
)


## ---- warning=FALSE, message=FALSE , fig.dim = c(6, 4)----
g.full <- buildSNNGraph(sce.pbmc, use.dimred = "PCA")
clust.full <- igraph::cluster_walktrap(g.full)$membership
sce.pbmc$clust.full <- factor(clust.full)
plotExpression(sce.pbmc,
    features = c("CD3E", "CCR7", "CD69", "CD44"),
    x = "clust.full", colour_by = "clust.full"
)


## ---- warning=FALSE, message=FALSE----------------------
# Repeating modelling and PCA on the subset of cells we have
# identified as memory T-cells (cluster 6).
memory <- 10
sce.memory <- sce.pbmc[, clust.full == memory]
dec.memory <- modelGeneVar(sce.memory)
sce.memory <- denoisePCA(sce.memory,
    technical = dec.memory,
    subset.row = getTopHVGs(dec.memory, prop = 0.1)
)
# Repeating clustering on the subset.
g.memory <- buildSNNGraph(sce.memory, use.dimred = "PCA")
clust.memory <- igraph::cluster_walktrap(g.memory)$membership
sce.memory$clust.memory <- factor(clust.memory)


## ---- warning=FALSE, message=FALSE, fig.dim = c(6, 4)----
plotExpression(sce.memory,
    features = c("CD8A", "CD4"),
    x = "clust.memory"
)


## -------------------------------------------------------
## Información de la sesión de R
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

