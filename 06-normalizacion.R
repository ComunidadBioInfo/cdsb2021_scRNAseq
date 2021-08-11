## ----echo=FALSE, fig.cap="Figura tomada de [1]"---------------------------------------
knitr::include_graphics("https://scrnaseq-course.cog.sanger.ac.uk/website/figures/RNA-Seq_workflow-5.pdf.jpg")


## ----message=FALSE, warning=FALSE-----------------------------------------------------
library("scRNAseq")
sce.zeisel <- ZeiselBrainData(ensembl = TRUE)
sce.zeisel


## ----message=FALSE, warning=FALSE, include=FALSE--------------------------------------
# Control de calidad
library("scater")

unfiltered <- sce.zeisel

is.mito <- which(rowData(sce.zeisel)$featureType == "mito")
stats <- perCellQCMetrics(sce.zeisel, subsets = list(Mt = is.mito))
qc <-
    quickPerCellQC(stats,
        percent_subsets = c("altexps_ERCC_percent", "subsets_Mt_percent")
    )

colSums(as.data.frame(qc))
sce.zeisel <- sce.zeisel[, !qc$discard]


## ----echo=FALSE-----------------------------------------------------------------------
knitr::include_graphics(here::here("img/libsize.png"))


## ----echo=FALSE-----------------------------------------------------------------------
knitr::include_graphics(here::here("img/libfactor.png"))


## -------------------------------------------------------------------------------------
# Estimar tamaños de librerías
lib.sf.zeisel <- librarySizeFactors(sce.zeisel)
# Examina la distribución de los tamaños de librerías
# que acabamos de estimar
summary(lib.sf.zeisel)
hist(log10(lib.sf.zeisel), xlab = "Log10[Library Size factor]", col = "grey80")


## ----echo=FALSE, fig.cap="Figura tomada de [2]", fig.height=5-------------------------
knitr::include_graphics("https://hbctraining.github.io/DGE_workshop/img/normalization_methods_composition.png")


## -------------------------------------------------------------------------------------
# Normalización por decircunvolución (deconvolution)
library("scran")
# Pre-clustering
set.seed(100)
clust.zeisel <- quickCluster(sce.zeisel)

# Calcula factores de tamaño para la decircunvolución (deconvolution)
deconv.sf.zeisel <-
    calculateSumFactors(sce.zeisel, clusters = clust.zeisel, min.mean = 0.1)
# Examina la distribución de los factores de tamaño
summary(deconv.sf.zeisel)
hist(log10(deconv.sf.zeisel),
    xlab = "Log10[Deconvolution size factor]",
    col = "grey80"
)


plot(lib.sf.zeisel,
    deconv.sf.zeisel,
    xlab = "Library size factor",
    ylab = "Deconvolution size factor",
    log = "xy",
    pch = 16
)
abline(a = 0, b = 1, col = "red")


## -------------------------------------------------------------------------------------
summary(deconv.sf.zeisel)


## -------------------------------------------------------------------------------------
50 - 10
1100 - 1000

log(50) - log(10)
log(1100) - log(1000)


## -------------------------------------------------------------------------------------
# Normalization
# set.seed(100)
# clust.zeisel <- quickCluster(sce.zeisel)
# sce.zeisel <- computeSumFactors(sce.zeisel, cluster=clust.zeisel, min.mean=0.1)

# Log transformation
sce.zeisel <- scater::logNormCounts(sce.zeisel)
assayNames(sce.zeisel)


## ---- message=FALSE, warning=FALSE----------------------------------------------------
library("Seurat")
# Create a Seurat obj
sce <- sce.zeisel
sce <- removeAltExps(sce)
seurat.zeisel <- as.Seurat(sce, counts = "counts", data = NULL)
seurat.zeisel

# Normalize using Seurat function
seurat.zeisel <- NormalizeData(seurat.zeisel, normalization.method = "LogNormalize")
# Compare Total counts per cell after normalization
ls.seurat <- colSums(seurat.zeisel[[SingleCellExperiment::mainExpName(x = sce)]]@data)
## Relacionado a
## https://github.com/satijalab/seurat/blob/9b3892961c9e1bf418af3bbb1bc79950adb481d7/R/objects.R#L1041-L1046
## donde podemos ver como Seurat convierte el objeto de SingleCellExperiment
## a un objeto de Seurat

summary(ls.seurat)
hist(ls.seurat)

# Trying to replicate it
ls.zeisel <- colSums(counts(sce.zeisel))
summary(ls.zeisel)
step1 <- t(counts(sce.zeisel)) / ls.zeisel # matrix(2,2,2,2) /c(1,2)
step2 <- step1 * 10000
step3 <- t(log1p(step2))
ls.steps <- colSums(step3)
summary(ls.steps)
plot(ls.seurat, ls.steps)

# Compare with deconv normalization
ls.log <- colSums(logcounts(sce.zeisel))
plot(ls.seurat, ls.log)


## ----echo=FALSE, fig.cap="Figura tomada de [2]", fig.height=5-------------------------
knitr::include_graphics("https://hbctraining.github.io/DGE_workshop/img/normalization_methods_depth.png")


## ----echo=FALSE, fig.cap="Figura tomada de [2]", fig.height=5-------------------------
knitr::include_graphics("https://hbctraining.github.io/DGE_workshop/img/normalization_methods_length.png")


## -------------------------------------------------------------------------------------
length(is.mito)


## -------------------------------------------------------------------------------------
colSums(as.data.frame(qc))


## -------------------------------------------------------------------------------------
# Plots
colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard

gridExtra::grid.arrange(
    plotColData(unfiltered, y = "sum", colour_by = "discard") +
        scale_y_log10() + ggtitle("Cuentas Totales"),
    plotColData(unfiltered, y = "detected", colour_by = "discard") +
        scale_y_log10() + ggtitle("Features (genes) detectados"),
    plotColData(unfiltered,
        y = "altexps_ERCC_percent",
        colour_by = "discard"
    ) + ggtitle("ERCC %"),
    plotColData(unfiltered,
        y = "subsets_Mt_percent",
        colour_by = "discard"
    ) + ggtitle("Mito %"),
    ncol = 2
)


## -------------------------------------------------------------------------------------
ls.zeisel <- colSums(counts(sce.zeisel))
summary(ls.zeisel)
hist(log10(ls.zeisel), xlab = "Log10[Library size]", col = "grey80")


## -------------------------------------------------------------------------------------
identical(lib.sf.zeisel, ls.zeisel)


## -------------------------------------------------------------------------------------
# Checamos proporcionalidad
plot(
    ls.zeisel,
    lib.sf.zeisel,
    log = "xy",
    main = "Proporcionalidad",
    xlab = "Library size",
    ylab = " Library size factor"
)


## -------------------------------------------------------------------------------------
## Ahora asegurate que su media sea 1 (unity mean)
lib_size_factors <- ls.zeisel / mean(ls.zeisel)
summary(lib_size_factors)
identical(lib_size_factors, lib.sf.zeisel)


## -------------------------------------------------------------------------------------
levels(clust.zeisel)


## -------------------------------------------------------------------------------------
cells_cluster <- sort(table(clust.zeisel))
cells_cluster
barplot(cells_cluster)


## -------------------------------------------------------------------------------------
set.seed(100)
sort(table(quickCluster(sce.zeisel, min.size = 200)))


## -------------------------------------------------------------------------------------
plot(lib.sf.zeisel,
    deconv.sf.zeisel,
    xlab = "Library size factor",
    ylab = "Deconvolution size factor",
    log = "xy",
    pch = 16,
    col = as.integer(factor(sce.zeisel$level1class))
)
abline(a = 0, b = 1, col = "red")
abline(a = -.2, b = 0.95, col = "red")
abline(a = 0.08, b = 1, col = "red")


## -------------------------------------------------------------------------------------
# sce.zeisel <- runPCA(sce.zeisel)
# plotPCA(sce.zeisel, colour_by = "level1class")
# plotRLE(sce.zeisel, exprs_values = "logcounts", colour_by = "level1class")


## -------------------------------------------------------------------------------------
## Información de la sesión de R
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

