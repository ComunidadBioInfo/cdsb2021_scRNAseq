## ---- warning=FALSE, message=FALSE--------------------------------------------------------
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


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
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


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
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


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
## Identificación de genes altamente variables
set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce.pbmc)
top.pbmc <- getTopHVGs(dec.pbmc, prop = 0.1)


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
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


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
# clustering
g <- buildSNNGraph(sce.pbmc, k = 10, use.dimred = "PCA")
clust <- igraph::cluster_walktrap(g)$membership
sce.pbmc$cluster <- factor(clust)


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
# Is gene 1 associated with the clustering?
plotExpression(sce.pbmc,
    features = rownames(sce.pbmc)[1],
    x = "cluster", colour_by = "cluster"
)


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
# Is gene 2 associated with the clustering?
plotExpression(sce.pbmc,
    features = rownames(sce.pbmc)[2],
    x = "cluster", colour_by = "cluster"
)


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
# Is gene 2512 associated with the clustering?
plotExpression(sce.pbmc,
    features = rownames(sce.pbmc)[2512],
    x = "cluster", colour_by = "cluster"
)


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
# Is gene CD3E associated with the clustering?
plotExpression(sce.pbmc,
    features = "CD3E",
    x = "cluster", colour_by = "cluster"
)


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
plotExpression(sce.pbmc,
    features = "CD3E",
    x = "cluster", colour_by = "cluster"
)


## ---- warning=FALSE, message=FALSE, echo=FALSE--------------------------------------------
library(kableExtra)
comparison <- c("1 vs 2", "1 vs 3", "...", "2 vs 1", "...", "18 vs 17")
logFC <- c("1.50", "-0.08", "...", "1.39", "...", "0.11")
Pval <- c("1.7e-198", "0.11", "...", "1.7e-198", "...", "0.46")
paired_tests <- data.frame(comparison, logFC, Pval)
knitr::kable(paired_tests, format = "html") %>%
    kable_styling(
        bootstrap_options = c("striped", "hover"),
        full_width = F,
        font_size = 12,
        position = "left"
    )


## ---- warning=FALSE, message=FALSE, eval = FALSE------------------------------------------
## scran::pairwiseTTests()
## scran::combineMarkers()


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
# scran::findMarkers()
library(scran)
markers.pbmc <- findMarkers(sce.pbmc,
    groups = sce.pbmc$cluster,
    test.type = "t", pval.type = "any"
)


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
chosen <- "9"
interesting <- markers.pbmc[[chosen]]


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
plotExpression(sce.pbmc, rownames(interesting)[1:4],
    x = "cluster", colour_by = "cluster"
)


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
best.set <- interesting[interesting$Top <= 6, ]
logFCs <- as.matrix(best.set[, -(1:3)])
colnames(logFCs) <- sub("logFC.", "", colnames(logFCs))
library(pheatmap)
pheatmap(logFCs, breaks = seq(-5, 5, length.out = 101))


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
markers.pbmc.up <- findMarkers(sce.pbmc,
    groups = sce.pbmc$cluster,
    test.type = "t", direction = "up", pval.type = "any"
)
interesting.up <- markers.pbmc.up[[chosen]]


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
markers.pbmc.up2 <- findMarkers(sce.pbmc,
    groups = sce.pbmc$cluster,
    test.type = "t", direction = "up", lfc = 1, pval.type = "any"
)
interesting.up2 <- markers.pbmc.up2[[chosen]]


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
best.set <- interesting.up2[interesting.up2$Top <= 5, ]
logFCs <- as.matrix(best.set[, -(1:3)])
colnames(logFCs) <- sub("logFC.", "", colnames(logFCs))
pheatmap(logFCs, breaks = seq(-5, 5, length.out = 101))


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
markers.pbmc.up3 <- findMarkers(sce.pbmc,
    groups = sce.pbmc$cluster,
    direction = "up", pval.type = "all"
)
interesting.up3 <- markers.pbmc.up3[[chosen]]


## ---- warning=FALSE, message=FALSE, echo=FALSE--------------------------------------------
library(kableExtra)
Poblacion <- c("DN(CD4-/CD8-)", "CD4+>", "CD8+>", "DP(CD4+/CD8+)")
Expresion_CD4 <- c("No", "Si", "No", "Si")
Expresion_CD8 <- c("No", "No", "Si", "Si")
poblacion <- data.frame(Poblacion, Expresion_CD4, Expresion_CD8)
knitr::kable(poblacion, format = "html") %>%
    kable_styling(
        bootstrap_options = c("striped", "hover"),
        full_width = F,
        font_size = 12,
        position = "left"
    )


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
markers.pbmc.up4 <- findMarkers(sce.pbmc,
    groups = sce.pbmc$cluster,
    direction = "up", pval.type = "some"
)
interesting.up4 <- markers.pbmc.up4[[chosen]]


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
markers.pbmc.wmw <- findMarkers(sce.pbmc,
    groups = sce.pbmc$cluster, test.type = "wilcox",
    direction = "up", pval.type = "any"
)
interesting.wmw <- markers.pbmc.wmw[[chosen]]


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
best.set <- interesting.wmw[interesting.wmw$Top <= 5, ]
AUCs <- as.matrix(best.set[, -(1:3)])
colnames(AUCs) <- sub("AUC.", "", colnames(AUCs))
pheatmap(AUCs,
    breaks = seq(0, 1, length.out = 21),
    color = viridis::viridis(21)
)


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
markers.pbmc.binom <- findMarkers(sce.pbmc,
    groups = sce.pbmc$cluster, test.type = "binom",
    direction = "up", pval.type = "any"
)
interesting.binom <- markers.pbmc.binom[[chosen]]


## ---- warning=FALSE, message=FALSE--------------------------------------------------------
top.genes <- head(rownames(interesting.binom))
plotExpression(sce.pbmc, x = "cluster", features = top.genes)


## -----------------------------------------------------------------------------------------
## Información de la sesión de R
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

