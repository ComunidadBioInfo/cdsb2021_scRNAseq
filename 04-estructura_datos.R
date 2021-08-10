## ---- warning=FALSE, message=FALSE------------------------------------------------------------------
library("scRNAseq")
sce.416b <- LunSpikeInData(which = "416b")

# Carga el paquete SingleCellExperiment
library("SingleCellExperiment")

"Primera parte
aquí checamos el slot assays"
# Extrae la matriz de cuentas del set de datos de 416b
counts.416b <- counts(sce.416b)
# CHEQUEMOS clase y dimensiones
class(counts.416b) # es matriz
dim(counts.416b) # indicará genes y células

# CONSTRUIR un nuevo objeto SCE de la matriz de cuentas !!!!!!
sce <- SingleCellExperiment(assays = list(counts = counts.416b))

# Revisa el objeto que acabamos de crear
sce

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB

# Accesa la matriz de cuenta del compartimento (slot) "assays"
# assays(sce, "counts")
# OJO: ¡esto puede inundar tu sesión de R!

# 1. El método general
assay(sce, "counts")[110:115, 1:3] # gene, cell
# 2. El método específico para accesar la matriz de cuentas "counts"
counts(sce)[110:115, 1:3]

# AGREGAR MAS ASSAYS
sce <- scater::logNormCounts(sce)
# Revisa el objeto que acabamos de actualizar
sce

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB

# 1. El método general
assay(sce, "logcounts")[110:115, 1:3]
# 2. El método específico para accesar la matriz de cuentas
# transformadas "logcounts"
logcounts(sce)[110:115, 1:3]


# agregemos un assay mas, esta vez de manera manual
assay(sce, "counts_100") <- assay(sce, "counts") + 100 # suma 100 a counts assay
# Enumera los "assays" en el objeto
assays(sce) # indica num y nombre de assays
assayNames(sce) # solo nos dará los nombres de los assays
#      assay(sce, "counts_100")[110:115, 1:3]
## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB


"segunda parte:
aquí checaremos metadata de las células"
# Extrae la información de las muestras (metadata) del set de datos de 416b
colData.416b <- colData(sce.416b) # podemos checar objeto en la cajita de environment de RStudio!!

# explorar datooos
table(colData.416b$phenotype)
table(colData.416b$block) # fue en varios dias?

# Agrega algo de esa información a nuestro objeto de SCE
colData(sce) <- colData.416b[, c("phenotype", "block")]
# Revisa el objeto que acabamos de actualizar
sce
# Accesa a la información de las muestras (metadata) en nuestro SCE
colData(sce) # usar head?
# Accesa una columna específica de la información de las muestras (metadata)
table(sce$block)
table(colData(sce)$block) # otra manera


# Ejemplo de una función que agrega columnas nuevas al colData
sce <- scater::addPerCellQC(sce.416b) # añade datos de control de calidad
# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
colData(sce)

# Revisa el objeto que acabamos de actualizar
sce

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB

## Agrega las cuentas normalizadas (lognorm) de nuevo
sce <- scater::logNormCounts(sce)

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB

# Ejemplo: obtén el subconjunto de células de fenotipo "wild type"
# Acuérdate que las células son columnas del SCE !!!!
sce[, sce$phenotype == "wild type phenotype"]



"Tercera parte:
examinaremos metadata de features (rowData)"
# Accesa la información de los genes de nuestro SCE
# ¡Está vació actualmente!
rowData(sce)

# Ejemplo de una función que agrega campos nuevos en el rowData
sce <- scater::addPerFeatureQC(sce)
# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
rowData(sce)

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB

# Descarga los archivos de anotación de la base de datos de Ensembl
# correspondientes usando los recursos disponibles vía AnnotationHub
library("AnnotationHub")
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))

# Obtén la posición del cromosoma para cada gen
ensdb <- ah[["AH73905"]]
chromosome <- mapIds(ensdb,
    keys = rownames(sce),
    keytype = "GENEID",
    column = "SEQNAME"
)
rowData(sce)$chromosome <- chromosome

# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
rowData(sce)

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB

# Ejemplo: obtén el subconjunto de datos donde los genes están en el
# cromosoma 3
# NOTA: which() fue necesario para lidear con los nombres de cromosoma
# que son NA
sce[which(rowData(sce)$chromosome == "3"), ]


"Cuarta parte:
examinamos slot metadata"
# Accesa la información de nuestro experimento usando metadata()
# ¡Está vació actualmente!
metadata(sce)

# La información en el metadata() es como Vegas - todo se vale
metadata(sce) <- list(
    favourite_genes = c("Shh", "Nck1", "Diablo"),
    analyst = c("Pete")
)

# Accesa la información de nuestro experimento usando metadata() de
# nuestro objeto actualizado
metadata(sce)

"Quinta parte:
examinamos slot de reducción de dimensiones"
# Ejemplo: agrega los componentes principales (PCs) de las logcounts
# NOTA: aprenderemos más sobre análisis de componentes principales (PCA) después
sce <- scater::runPCA(sce)
# Revisa el objeto que acabamos de actualizar
sce
# Accesa la matriz de PCA del componente (slot) reducedDims
reducedDim(sce, "PCA")[1:6, 1:3]

# Ejemplo, agrega una representación de los logcounts en t-SNE
# NOTA: aprenderemos más sobre t-SNE después
sce <- scater::runTSNE(sce)
# Revisa el objeto que acabamos de actualizar
sce
# Accesa a la matriz de t-SNE en el componente (slot) de reducedDims
head(reducedDim(sce, "TSNE"))

# Ejemplo: agrega una representación 'manual' de los logcounts en UMAP
# NOTA: aprenderemos más sobre UMAP después y de una forma más sencilla de
#       calcularla
u <- uwot::umap(t(logcounts(sce)), n_components = 2)

# Agrega la matriz de UMAP al componente (slot) reducedDims
reducedDim(sce, "UMAP") <- u

# Accesa a la matriz de UMAP desde el componente (slot) reducedDims
head(reducedDim(sce, "UMAP"))

# Enumera los resultados de reducción de dimensiones en nuestro objeto SCE
reducedDims(sce)

"Sexta parte:
experimentos alternativos"
# Extrae la información de ERCC de nuestro SCE para el set de datos de 416b
ercc.sce.416b <- altExp(sce.416b, "ERCC")
# Inspecciona el SCE para los datos de ERCC
ercc.sce.416b

# Agrega el SCE de ERCC como un experimento alternativo a nuestro SCE
altExp(sce, "ERCC") <- ercc.sce.416b
# Revisa el objeto que acabamos de actualizar
sce

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce) / 1024^2 ## En MB

# Enumera los experimentos alternativos almacenados en nuestro objeto
altExps(sce)


# El crear un subconjunto del SCE por muestra (célula) automáticamente
# obtiene el subconjunto de los experimentos alternativos
sce.subset <- sce[, 1:10]
ncol(sce.subset)
ncol(altExp(sce.subset))

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce.subset) / 1024^2 ## En MB

"Septima parte:
sizefactores"
# Extrae los factores de tamaño (size factors)
# Estos fueron añadidos a nuestro objeto cuando corrimos
# scater::logNormCounts(sce)
head(sizeFactors(sce))

# "Automáticamente" reemplaza los factores de tamaño
sce <- scran::computeSumFactors(sce)
head(sizeFactors(sce))

# "Manualmente" reemplaza los factores de tamaño
sizeFactors(sce) <- scater::librarySizeFactors(sce)
head(sizeFactors(sce))


## ---------------------------------------------------------------------------------------------------
# Creamos un data.frame
df <- data.frame(x = c(TRUE, FALSE, NA, NA), y = c(12, 34, 56, 78))
row.names(df) <- letters[1:4]

df

# Para acceder a los nombres de las columnas
colnames(df)

# Para acceder a los nombres de las filas
rownames(df)

# Podemos sacar información booleana
df$y < 20

# Y podemos acceder al mismo data.frame con la información booleana

df[df$y < 40, ]

## %in% (dentro de)
bool_info <- rownames(df) %in% c("a", "c", "z")
df[bool_info, ]

## & (y)
bool_info <- df$y < 50 & df$y > 20
df[bool_info, ]

## | (o)
bool_info <- df$y < 20 | df$y > 60
df[bool_info, ]


## ---------------------------------------------------------------------------------------------------
library("SingleCellExperiment")
library("scRNAseq")

# Mini muestreo del set de datos usado en: https://bioconductor.org/books/release/OSCA/zeisel-mouse-brain-strt-seq.html#introduction-5

archivo_cuentas <- "https://raw.githubusercontent.com/emarquezz/minidataset_osca/main/min_sce.csv"
archivo_rowData <- "https://raw.githubusercontent.com/emarquezz/minidataset_osca/main/rowD.csv"
archivo_colData <- "https://raw.githubusercontent.com/emarquezz/minidataset_osca/main/colD.csv"


counts <- read.csv(archivo_cuentas, row.names = 1, header = TRUE, check.names = F)
col.data <- DataFrame(read.csv(archivo_colData, row.names = 1, header = TRUE, check.names = F))
row.data <- read.csv(archivo_rowData, row.names = 1, header = TRUE, check.names = F)


## ---------------------------------------------------------------------------------------------------
int_gen <- c("Angpt1", "Chic2", "Mir503", "Magee2", "Nenf", "Eps15l1", "Hsf2bp", "Gnptg", "Vegfb", "Atmin", "Gad1", "Gad2", "Slc32a1", "Dner", "Slc2a13", "Slc6a1", "Nrxn3")


## ---- eval= FALSE-----------------------------------------------------------------------------------
## library("scater")
## 
## plotHeatmap(object = tej_min_sce, features = rownames(tej_min_sce), order_columns_by = "level1class")


## ---------------------------------------------------------------------------------------------------
sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = col.data,
    rowData = row.data
)

sce <- scater::logNormCounts(sce)

sce


## ---------------------------------------------------------------------------------------------------
bool_data <- rownames(rowData(sce))
min_sce <- sce[bool_data %in% int_gen, ]


## ---------------------------------------------------------------------------------------------------
tej_int <- min_sce$level1class == "interneurons"
tej_pyr <- min_sce$level1class == "pyramidal CA1"


tej_min_sce <- min_sce[, tej_int | tej_pyr]


## ---------------------------------------------------------------------------------------------------
library("scater")

plotHeatmap(object = tej_min_sce, features = rownames(tej_min_sce), order_columns_by = "level1class")


## ---------------------------------------------------------------------------------------------------
int_gen <- c("ENSMUSG00000071076", "ENSMUSG00000002076", "ENSMUSG00000024962", "ENSMUSG00000031224", "ENSMSG00000036560", "ENSMUSG00000037499", "ENSMUSG00000006276", "ENSMUSG00000035521", "ENSMUSG00000047388", "ENSMUSG0000051079", "ENSMUSG00000076122", "ENSMUSG00000029229", "ENSMUSG00000022309", "ENSMUSG00000036766", "ENSMUSG00000070880", "ENSMUSG00000026787", "ENSMUSG00000066392", "ENSMUSG00000036298", "ENSMUSG00000037771", "ENSMUSG00000030310")


## ---------------------------------------------------------------------------------------------------
# Descarga datos de ejemplo procesados con CellRanger
# Paréntesis: al usar BiocFileCache solo tenemos que descargar
#             los datos una vez.
library("BiocFileCache")
bfc <- BiocFileCache()
pbmc.url <-
    paste0(
        "http://cf.10xgenomics.com/samples/cell-vdj/",
        "3.1.0/vdj_v1_hs_pbmc3/",
        "vdj_v1_hs_pbmc3_filtered_feature_bc_matrix.tar.gz"
    )
pbmc.data <- bfcrpath(bfc, pbmc.url)

# Extrae los archivos en un directorio temporal
untar(pbmc.data, exdir = tempdir())

# Enumera los archivos que descargamos y que extrajimos
# Estos son los archivos típicos de CellRanger
pbmc.dir <- file.path(
    tempdir(),
    "filtered_feature_bc_matrix"
)
list.files(pbmc.dir)

# Importa los datos como un objeto de tipo SingleCellExperiment
library("DropletUtils")
sce.pbmc <- read10xCounts(pbmc.dir)
# Revisa el objeto que acabamos de construir
sce.pbmc

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce.pbmc) / 1024^2 ## En MB

# Almacena la información de CITE-seq como un experimento alternativo
sce.pbmc <- splitAltExps(sce.pbmc, rowData(sce.pbmc)$Type)
# Revisa el objeto que acabamos de actualizar
sce.pbmc

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce.pbmc) / 1024^2 ## En MB

# Descarga datos de ejemplo procesados con scPipe
library("BiocFileCache")
bfc <- BiocFileCache()
sis_seq.url <-
    "https://github.com/LuyiTian/SIS-seq_script/archive/master.zip"
sis_seq.data <- bfcrpath(bfc, sis_seq.url)

# Extrae los archivos en un directorio temporal
unzip(sis_seq.data, exdir = tempdir())

# Enumera (algunos de) los archivos que descargamos y extrajimos
# Estos son los archivos típicos de scPipe
sis_seq.dir <- file.path(
    tempdir(),
    "SIS-seq_script-master",
    "data",
    "BcorKO_scRNAseq",
    "RPI10"
)
list.files(sis_seq.dir)

# Importa los datos como un objeto de tipo SingleCellExperiment
library("scPipe")
sce.sis_seq <- create_sce_by_dir(sis_seq.dir)
# Revisa el objeto que acabamos de construir
sce.sis_seq

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(sce.sis_seq) / 1024^2 ## En MB

# Descarga un ejemplo de un montón de archivos
library("BiocFileCache")
bfc <- BiocFileCache()
lun_counts.url <-
    paste0(
        "https://www.ebi.ac.uk/arrayexpress/files/",
        "E-MTAB-5522/E-MTAB-5522.processed.1.zip"
    )
lun_counts.data <- bfcrpath(bfc, lun_counts.url)
lun_coldata.url <-
    paste0(
        "https://www.ebi.ac.uk/arrayexpress/files/",
        "E-MTAB-5522/E-MTAB-5522.sdrf.txt"
    )
lun_coldata.data <- bfcrpath(bfc, lun_coldata.url)

# Extrae los archivos en un directorio temporal
lun_counts.dir <- tempfile("lun_counts.")
unzip(lun_counts.data, exdir = lun_counts.dir)

# Enumera los archivos que descargamos y extrajimos
list.files(lun_counts.dir)

# Lee la matriz de cuentas (para una placa)
lun.counts <- read.delim(
    file.path(lun_counts.dir, "counts_Calero_20160113.tsv"),
    header = TRUE,
    row.names = 1,
    check.names = FALSE
)
# Almacena la información de la longitud de los genes para después
gene.lengths <- lun.counts$Length
# Convierte los datos de cuentas de genez a una matriz (quitamos las longitudes)
lun.counts <- as.matrix(lun.counts[, -1])

# Lee la información de las muestras (células)
lun.coldata <- read.delim(lun_coldata.data,
    check.names = FALSE,
    stringsAsFactors = FALSE
)
library("S4Vectors")
lun.coldata <- as(lun.coldata, "DataFrame")

# Pon en orden la información de las muestras para que
# sea idéntico al orden en la matriz de cuentas
m <- match(
    colnames(lun.counts),
    lun.coldata$`Source Name`
)
lun.coldata <- lun.coldata[m, ]

# Construye la tabla de información de los genes
lun.rowdata <- DataFrame(Length = gene.lengths)

# Construye el objeto de SingleCellExperiment
lun.sce <- SingleCellExperiment(
    assays = list(assays = lun.counts),
    colData = lun.coldata,
    rowData = lun.rowdata
)
# Revisa el objeto que acabamos de construir
lun.sce

## ¿Qué tan grande es el objeto de R?
lobstr::obj_size(lun.sce) / 1024^2 ## En MB


## ---------------------------------------------------------------------------------------------------
## Información de la sesión de R
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

