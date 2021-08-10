## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
## Cargar paquetes de R
library("BiocFileCache") ## para descargar datos
library("dplyr") ## para filtar datos
library("Seurat") ## paquete principal de este capítulo
library("patchwork") ## para graficar imágenes juntas


## ----------------------------------------------------------------------------------------------------
# Usemos datos de pbmc3k tal y como lo hacen en
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
# pero con nuestro propio código
bfc <- BiocFileCache()
raw.path <- bfcrpath(bfc, file.path(
    "http://cf.10xgenomics.com/samples",
    "cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
))
untar(raw.path, exdir = file.path(tempdir(), "pbmc3k"))
fname <- file.path(tempdir(), "pbmc3k/filtered_gene_bc_matrices/hg19")

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = fname)
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
str(pbmc)


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size
dense.size / sparse.size


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
# Filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


## ---- warning=FALSE, message=FALSE-------------------------------------------------------------------
head(pbmc@meta.data, 5)


## ----------------------------------------------------------------------------------------------------
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


## ---- fig.width = 14, fig.height = 7-----------------------------------------------------------------
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


## ----------------------------------------------------------------------------------------------------
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


## ----------------------------------------------------------------------------------------------------
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))


## ----------------------------------------------------------------------------------------------------
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")


## ----------------------------------------------------------------------------------------------------
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


## ----------------------------------------------------------------------------------------------------
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)


## ----------------------------------------------------------------------------------------------------
JackStrawPlot(pbmc, dims = 1:15)


## ----------------------------------------------------------------------------------------------------
ElbowPlot(pbmc)


## ----------------------------------------------------------------------------------------------------
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)


## ----------------------------------------------------------------------------------------------------
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")


## ----------------------------------------------------------------------------------------------------
if (interactive()) {
    saveRDS(pbmc, file = "pbmc_tutorial.rds")
}


## ----------------------------------------------------------------------------------------------------
# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_log2FC)


## ----------------------------------------------------------------------------------------------------
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


## ----------------------------------------------------------------------------------------------------
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

## you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c(
    "MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"
))

# DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


## ----------------------------------------------------------------------------------------------------
new.cluster.ids <- c(
    "Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet"
)
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


## ----------------------------------------------------------------------------------------------------
if (interactive()) {
    saveRDS(pbmc, file = "pbmc3k_final.rds")
}


## ----------------------------------------------------------------------------------------------------
## Información de la sesión de R
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

