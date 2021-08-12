## ----cargar_paquetes, message = FALSE-------------------
## Paquetes de este capítulo
library("MouseGastrulationData") ## para descargar datos de ejemplo
library("scater") ## para gráficas y control de calidad
library("scran") ## para selección de genes, clustering, etc
library("batchelor") ## para métodos de correción de batch (lote)
library("patchwork") ## para agrupar gráficas
library("Polychrome") ## para muchos colores
library("bluster") ## métodos de clustering
library("edgeR") ## para expresión diferencial


## ---- echo = FALSE--------------------------------------
set.seed(20210811)
knitr::kable(data.frame(
    gene = paste0("gene", rep(rep(1:2, each = 2), 2)),
    condición = paste0("grupo", rep(rep(1:2, 2), 2)),
    celula = paste0("celula", rep(1:2, each = 4)),
    expresión = round(rnorm(8, mean = 10, sd = 2), 2)
))


## ---- echo = FALSE--------------------------------------
set.seed(20210811)
knitr::kable(data.frame(
    condición = paste0("grupo", rep(1:2, 2)),
    celula = paste0("celula", rep(1:2, each = 2)),
    frecuencia = round(rnorm(4, mean = 40, sd = 4), 0)
))


## -------------------------------------------------------
#--- loading ---#
library("MouseGastrulationData")
sce.chimera <- WTChimeraData(samples = 5:10)
sce.chimera


## -------------------------------------------------------
## Exploremos los datos
sapply(colData(sce.chimera)[, -(1:2)], function(x) {
    x <- if (is.character(x) || is.integer(x)) factor(x) else x
    summary(x)
})


## -------------------------------------------------------
## Número de células en nuestras variables principales
with(colData(sce.chimera), table(sample, pool, tomato))

## Número de tipos celulares
length(unique(sce.chimera$celltype.mapped))


## -------------------------------------------------------
#--- feature-annotation ---#
library("scater")
rownames(sce.chimera) <- uniquifyFeatureNames(
    rowData(sce.chimera)$ENSEMBL, rowData(sce.chimera)$SYMBOL
)

#--- quality-control ---#
drop <- sce.chimera$celltype.mapped %in% c("stripped", "Doublet")
sce.chimera <- sce.chimera[, !drop]

#--- normalization ---#
sce.chimera <- logNormCounts(sce.chimera)

#--- variance-modelling ---#
library("scran")
dec.chimera <- modelGeneVar(sce.chimera, block = sce.chimera$sample)
chosen.hvgs <- dec.chimera$bio > 0

#--- merging ---#
library("batchelor")
set.seed(01001001)
merged <- correctExperiments(sce.chimera,
    batch = sce.chimera$sample,
    subset.row = chosen.hvgs,
    PARAM = FastMnnParam(
        merge.order = list(
            list(1, 3, 5), # WT (3 replicates)
            list(2, 4, 6) # td-Tomato (3 replicates)
        )
    )
)

#--- clustering ---#
g <- buildSNNGraph(merged, use.dimred = "corrected")
clusters <- igraph::cluster_louvain(g)
colLabels(merged) <- factor(clusters$membership)

#--- dimensionality-reduction ---#
merged <- runTSNE(merged, dimred = "corrected", external_neighbors = TRUE)
merged <- runUMAP(merged, dimred = "corrected", external_neighbors = TRUE)


## -------------------------------------------------------
## Clusters vs DE por td-Tomato
table(colLabels(merged), merged$tomato)

## Clusters vs lotes de muestras (batch)
table(colLabels(merged), merged$pool)


## ---- fig.width=10--------------------------------------
library("patchwork")
plotTSNE(merged, colour_by = "tomato", text_by = "label") +
    plotTSNE(merged, colour_by = data.frame(pool = factor(merged$pool)))


## ---- fig.width=10, warning = FALSE---------------------
## Definir colores, si no scater nos los pone en una escala
## continua
cols_label <- Polychrome::palette36.colors(length(unique(merged$label)))
names(cols_label) <- unique(merged$label)
cols_celltype.mapped <- Polychrome::palette36.colors(length(unique(merged$celltype.mapped)))
names(cols_celltype.mapped) <- unique(merged$celltype.mapped)

## Nuestros clusters vs anotación de células por los
## autores originales
plotTSNE(merged, colour_by = "label", text_by = "label") +
    theme(legend.position = "none") +
    scale_colour_manual(values = cols_label) +
    plotTSNE(merged, colour_by = "celltype.mapped") +
    theme(legend.position = "none") +
    scale_colour_manual(values = cols_celltype.mapped)


## -------------------------------------------------------
library("bluster")
pairwiseRand(colLabels(merged), merged$celltype.mapped, "index")

by.label <- table(colLabels(merged), merged$celltype.mapped)
pheatmap::pheatmap(log2(by.label + 1), color = viridis::viridis(101))


## -------------------------------------------------------
# Using 'label' and 'sample' as our two factors; each column of the output
# corresponds to one unique combination of these two factors.
summed <- aggregateAcrossCells(merged,
    id = colData(merged)[, c("celltype.mapped", "sample")]
)
summed

dim(merged)
dim(summed)

with(colData(merged), length(unique(celltype.mapped)) * length(unique(sample)))


## -------------------------------------------------------
label <- "Mesenchyme"
current <- summed[, label == summed$celltype.mapped]
dim(current)


## -------------------------------------------------------
# Creating up a DGEList object for use in edgeR:
library("edgeR")
y <- DGEList(counts(current), samples = colData(current))
y


## -------------------------------------------------------
discarded <- current$ncells < 10
y <- y[, !discarded]
summary(discarded)


## -------------------------------------------------------
keep <- filterByExpr(y, group = current$tomato)
y <- y[keep, ]
summary(keep)


## -------------------------------------------------------
y <- calcNormFactors(y)
y$samples


## -------------------------------------------------------
par(mfrow = c(2, 3))
for (i in seq_len(ncol(y))) {
    plotMD(y, column = i)
}


## -------------------------------------------------------
par(mfrow = c(1, 1))
plotMDS(cpm(y, log = TRUE),
    col = ifelse(y$samples$tomato, "red", "blue")
)


## -------------------------------------------------------
design <- model.matrix(~ factor(pool) + factor(tomato), y$samples)
design


## -------------------------------------------------------
if (interactive()) {
    ExploreModelMatrix::ExploreModelMatrix(y$samples[, c("pool", "tomato")], ~ factor(pool) + factor(tomato))
}


## -------------------------------------------------------
y <- estimateDisp(y, design)
summary(y$trended.dispersion)

plotBCV(y)

fit <- glmQLFit(y, design, robust = TRUE)
summary(fit$var.prior)
summary(fit$df.prior)
plotQLDisp(fit)


## -------------------------------------------------------
res <- glmQLFTest(fit, coef = ncol(design))
summary(decideTests(res))

topTags(res)


## -------------------------------------------------------
# Removing all pseudo-bulk samples with 'insufficient' cells.
summed.filt <- summed[, summed$ncells >= 10]

library("scran")
de.results <- pseudoBulkDGE(summed.filt,
    label = summed.filt$celltype.mapped,
    design = ~ factor(pool) + tomato,
    coef = "tomatoTRUE",
    condition = summed.filt$tomato
)
class(de.results)
length(de.results)


## -------------------------------------------------------
cur.results <- de.results[["Allantois"]]
cur.results[order(cur.results$PValue), ]

y.allantois <- metadata(cur.results)$y
plotBCV(y.allantois)


## -------------------------------------------------------
metadata(de.results)$failed


## -------------------------------------------------------
cur.results.Mesenchyme <- de.results[["Mesenchyme"]]
y.Mesenchyme <- metadata(cur.results.Mesenchyme)$y
plotBCV(y.Mesenchyme)


## -------------------------------------------------------
table(merged$sample)


## -------------------------------------------------------
abundances <- table(merged$celltype.mapped, merged$sample)
abundances <- unclass(abundances)
head(abundances)


## -------------------------------------------------------
# Attaching some column metadata.
extra.info <- colData(merged)[match(colnames(abundances), merged$sample), ]
y.ab <- DGEList(abundances, samples = extra.info)
y.ab


## -------------------------------------------------------
keep <- filterByExpr(y.ab, group = y.ab$samples$tomato)
y.ab <- y.ab[keep, ]
summary(keep)


## -------------------------------------------------------
design <- model.matrix(~ factor(pool) + factor(tomato), y.ab$samples)
design

y.ab <- estimateDisp(y.ab, design, trend = "none")
summary(y.ab$common.dispersion)

plotBCV(y.ab, cex = 1)


## -------------------------------------------------------
fit.ab <- glmQLFit(y.ab, design, robust = TRUE, abundance.trend = FALSE)
summary(fit.ab$var.prior)

plotQLDisp(fit.ab, cex = 1)


## -------------------------------------------------------
res <- glmQLFTest(fit.ab, coef = ncol(design))
summary(decideTests(res))

topTags(res)


## -------------------------------------------------------
## Información de la sesión de R
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

