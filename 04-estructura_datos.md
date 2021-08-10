# Estructura e importe de datos

Instructoras: [Elisa Márquez Zavala](https://twitter.com/naielisha), [Citlali Gil Aguillon](http://twitter.com/argininaa)

Contenido adaptado de [CDSB2020: Introducción a scRNA-seq, estructura e importe de datos
](https://comunidadbioinfo.github.io/cdsb2020/scRNAseq/02-data-infrastructure-and-import.html#1) de [Leonardo Collado Torres](https://github.com/lcolladotor) 


## Estructura de `SingleCellExperiment`

Dentro del curso -y más ampliamente dentro de los análisis de scRNA-seq en R- emplearemos la clase `SingleCellExperiment` . 

>Además recordemos que la interoperabilidad dentro de los paquetes de Bioconductor hará que facilmente puedas ajustarte más facilmente a otros paquetes que vayamos encontrando útiles (la infraestructura de los dato nos seguirá sirviendo!).

Podríamos dividir esta clase en cuatro categorias:

- datos primarios y transformados (donde estaran lo)

- metadata de datos (información de los genes o features, de las células y del experimento)

- reduccion de dimensiones 

- experimentos alternativos

Todo lo anterior está ligado, por lo que hace más sencillo el manejo de subsets de interés

> por ejemplo si nos interesa los genes x,y,z podríamos solicitarlos y con ello traer la información de las demás tablas

Examinaremos cada una de estas partes a detalle


Usaremos las diapositivas de [Peter Hickey](https://www.peterhickey.org/) para explicar la clase `SingleCellExperiment` [LINK A DIAPOSITIVAS](https://docs.google.com/presentation/d/1X9qP3wNlnn3BMUQhuZwAo4vCV76c33X_M-UnHxkPZpE/edit#slide=id.g7e63e0fe24_0_151). A la par se demostrará el siguiente código


## Ejercicio 1

The A. T. L. Lun et al. (2017) dataset contains two 96-well plates of 416B cells (an immortalized mouse myeloid progenitor cell line), processed using the Smart-seq2 protocol (Picelli et al. 2014). A constant amount of spike-in RNA from the External RNA Controls Consortium (ERCC) was also added to each cell’s lysate prior to library preparation. High-throughput sequencing was performed and the expression of each gene was quantified by counting the total number of reads mapped to its exonic regions. Similarly, the quantity of each spike-in transcript was measured by counting the number of reads mapped to the spike-in reference sequences.

Fragmento obtenido de https://bioconductor.org/books/release/OSCA/lun-416b-cell-line-smart-seq2.html


```r
library('scRNAseq')
sce.416b <- LunSpikeInData(which = "416b")

# Carga el paquete SingleCellExperiment
library('SingleCellExperiment')

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
pryr::object_size(sce)

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
pryr::object_size(sce)

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
pryr::object_size(sce)


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
table(colData(sce)$block) #otra manera


# Ejemplo de una función que agrega columnas nuevas al colData
sce <- scater::addPerCellQC(sce.416b) # añade datos de control de calidad
# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
colData(sce)

# Revisa el objeto que acabamos de actualizar
sce

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce)

## Agrega las cuentas normalizadas (lognorm) de nuevo
sce <- scater::logNormCounts(sce)

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce)

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
pryr::object_size(sce)

# Descarga los archivos de anotación de la base de datos de Ensembl
# correspondientes usando los recursos disponibles vía AnnotationHub
library('AnnotationHub')
ah <- AnnotationHub()
query(ah, c("Mus musculus", "Ensembl", "v97"))

# Obtén la posición del cromosoma para cada gen
ensdb <- ah[["AH73905"]]
chromosome <- mapIds(ensdb,
    keys = rownames(sce),
    keytype = "GENEID",
    column = "SEQNAME")
rowData(sce)$chromosome <- chromosome

# Accesa a la información de las muestras (metadata) en nuestro SCE actualizado
rowData(sce)

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce)

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
metadata(sce) <- list(favourite_genes = c("Shh", "Nck1", "Diablo"),
    analyst = c("Pete"))

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
pryr::object_size(sce)

# Enumera los experimentos alternativos almacenados en nuestro objeto
altExps(sce)


# El crear un subconjunto del SCE por muestra (célula) automáticamente
# obtiene el subconjunto de los experimentos alternativos
sce.subset <- sce[, 1:10]
ncol(sce.subset)
ncol(altExp(sce.subset))

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce.subset)

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
```

NOTA: La clase `SingleCellExperiment`está basada en `SummarizedExperiment`, por lo que ya estamos un poco familiarizados con esta nueva clase.

> Una de las diferencias es que no contiene ranuras para análisis de reducción de dimensiones. (será solo datos y metadata de ellos)

### Preguntas:

* ¿Cuáles son los tipos de tablas que debe siempre contenter el objeto `sce`?

* ¿Donde usamos los `colnames(sce)`?

* ¿donde usamos los `rownames(sce)`?

## Breve repaso de R

```r
# Creamos un data.frame
df <- data.frame(x = c(TRUE, FALSE, NA, NA), y = c(12, 34, 56, 78))
row.names(df) <- letters[1:4]

df 
```

```
##       x  y
## a  TRUE 12
## b FALSE 34
## c    NA 56
## d    NA 78
```

```r
# Para acceder a los nombres de las columnas
colnames(df)
```

```
## [1] "x" "y"
```

```r
# Para acceder a los nombres de las filas
rownames(df)
```

```
## [1] "a" "b" "c" "d"
```

```r
# Podemos sacar información booleana
df$y < 20
```

```
## [1]  TRUE FALSE FALSE FALSE
```

```r
# Y podemos acceder al mismo data.frame con la información booleana

df[df$y<40,]
```

```
##       x  y
## a  TRUE 12
## b FALSE 34
```

```r
## %in% (dentro de)
bool_info <- rownames(df) %in% c('a','c','z')
df[bool_info,]
```

```
##      x  y
## a TRUE 12
## c   NA 56
```

```r
## & (y)  
bool_info <- df$y < 50 & df$y > 20
df[bool_info,]
```

```
##       x  y
## b FALSE 34
```

```r
## | (o)
bool_info <- df$y < 20 | df$y > 60
df[bool_info,]
```

```
##      x  y
## a TRUE 12
## d   NA 78
```


## Ejercicio 2


```r
library("SingleCellExperiment")
library("scRNAseq")

# Mini muestreo del set de datos usado en: https://bioconductor.org/books/release/OSCA/zeisel-mouse-brain-strt-seq.html#introduction-5

archivo_cuentas <- "https://raw.githubusercontent.com/emarquezz/minidataset_osca/main/min_sce.csv"
archivo_rowData <-"https://raw.githubusercontent.com/emarquezz/minidataset_osca/main/rowD.csv"
archivo_colData <- "https://raw.githubusercontent.com/emarquezz/minidataset_osca/main/colD.csv"


counts <- read.csv(archivo_cuentas, row.names = 1, header= TRUE, check.names=F)
col.data <- DataFrame(read.csv(archivo_colData, row.names = 1, header= TRUE, check.names=F)) 
row.data <- read.csv(archivo_rowData, row.names = 1, header= TRUE, check.names=F)
```


* Crea un objeto SingleCellExperiment


  - ¿Cuántos genes tenemos?

  - ¿Qué información tenemos en el rowData? 

* Extraigan los datos de los genes que nos interesan (objeto ´int_gen´)
  - Pasen el mouse sobre los siguientes textos para ver recomendaciones:
  - [Recomendación 1](## "Recuerda que puedes acceder con booleanos")
  
  
  - [Recomendación 2](## "Recuerda que puedes utilizar el operador %in%")
  
  
  - [Recomendación 3](## "Esto es como el ejercicio 'Using **rowData** for subsetting'")


```r
int_gen <-  c("Angpt1","Chic2","Mir503","Magee2","Nenf","Eps15l1","Hsf2bp","Gnptg","Vegfb","Atmin","Gad1","Gad2","Slc32a1","Dner","Slc2a13","Slc6a1","Nrxn3")
```

* Creen un objeto llamado ´min_sce´ con los datos de los genes



  - ¿Cuáles son parte del tejido ´interneurons´ o del tejido ´pyramidal CA1´ ? (del objeto ´min_sce´) 


  - [Recomendación 1](## "Recuerda que puedes acceder con booleanos")
  
  
  - [Recomendación 2](## "Recuerda que puedes utilizar el operador |")
  
  
  - [Recomendación 3](## "Esto es como el ejercicio 'Using **colData** for subsetting'")


* Con este subconjunto, crea el objeto tej_min_sce



* Una vez que tengan el objeto ´SingleCellExperiment´ llamado ´tej_min_sce´, corran el siguiente código.

```r
library("scater")

plotHeatmap(object= tej_min_sce, features=rownames(tej_min_sce), order_columns_by="level1class")
```

### Extra

* Realiza los mismos pasos, pero ahora los genes que buscamos no tienen el nombre usual de gen (Gad1), sino su Ensembl gene IDs


```r
int_gen <-  c("ENSMUSG00000071076","ENSMUSG00000002076","ENSMUSG00000024962","ENSMUSG00000031224","ENSMSG00000036560","ENSMUSG00000037499","ENSMUSG00000006276","ENSMUSG00000035521","ENSMUSG00000047388","ENSMUSG0000051079","ENSMUSG00000076122","ENSMUSG00000029229","ENSMUSG00000022309","ENSMUSG00000036766","ENSMUSG00000070880","ENSMUSG00000026787","ENSMUSG00000066392","ENSMUSG00000036298","ENSMUSG00000037771","ENSMUSG00000030310")
```

## Import data

* Diapositivas de Peter Hickey

Ve las diapositivas [aquí](https://docs.google.com/presentation/d/1X9qP3wNlnn3BMUQhuZwAo4vCV76c33X_M-UnHxkPZpE/edit)


```r
# Descarga datos de ejemplo procesados con CellRanger
# Paréntesis: al usar BiocFileCache solo tenemos que descargar
#             los datos una vez.
library('BiocFileCache')
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
pbmc.dir <- file.path(tempdir(),
    "filtered_feature_bc_matrix")
list.files(pbmc.dir)

# Importa los datos como un objeto de tipo SingleCellExperiment
library('DropletUtils')
sce.pbmc <- read10xCounts(pbmc.dir)
# Revisa el objeto que acabamos de construir
sce.pbmc

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce.pbmc)

# Almacena la información de CITE-seq como un experimento alternativo
sce.pbmc <- splitAltExps(sce.pbmc, rowData(sce.pbmc)$Type)
# Revisa el objeto que acabamos de actualizar
sce.pbmc

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce.pbmc)

# Descarga datos de ejemplo procesados con scPipe
library('BiocFileCache')
bfc <- BiocFileCache()
sis_seq.url <-
    "https://github.com/LuyiTian/SIS-seq_script/archive/master.zip"
sis_seq.data <- bfcrpath(bfc, sis_seq.url)

# Extrae los archivos en un directorio temporal
unzip(sis_seq.data, exdir = tempdir())

# Enumera (algunos de) los archivos que descargamos y extrajimos
# Estos son los archivos típicos de scPipe
sis_seq.dir <- file.path(tempdir(),
    "SIS-seq_script-master",
    "data",
    "BcorKO_scRNAseq",
    "RPI10")
list.files(sis_seq.dir)

# Importa los datos como un objeto de tipo SingleCellExperiment
library('scPipe')
sce.sis_seq <- create_sce_by_dir(sis_seq.dir)
# Revisa el objeto que acabamos de construir
sce.sis_seq

## ¿Qué tan grande es el objeto de R?
pryr::object_size(sce.sis_seq)

# Descarga un ejemplo de un montón de archivos
library('BiocFileCache')
bfc <- BiocFileCache()
lun_counts.url <-
    paste0(
        "https://www.ebi.ac.uk/arrayexpress/files/",
        "E-MTAB-5522/E-MTAB-5522.processed.1.zip"
    )
lun_counts.data <- bfcrpath(bfc, lun_counts.url)
lun_coldata.url <-
    paste0("https://www.ebi.ac.uk/arrayexpress/files/",
        "E-MTAB-5522/E-MTAB-5522.sdrf.txt")
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
    stringsAsFactors = FALSE)
library('S4Vectors')
lun.coldata <- as(lun.coldata, "DataFrame")

# Pon en orden la información de las muestras para que
# sea idéntico al orden en la matriz de cuentas
m <- match(colnames(lun.counts),
    lun.coldata$`Source Name`)
lun.coldata <- lun.coldata[m,]

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
pryr::object_size(lun.sce)
```

## Detalles de la sesión de R


```r
## Información de la sesión de R
Sys.time()
```

```
## [1] "2021-08-10 09:11:30 UTC"
```

```r
proc.time()
```

```
##    user  system elapsed 
##   0.552   0.128   0.560
```

```r
options(width = 120)
sessioninfo::session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
##  setting  value                       
##  version  R version 4.1.0 (2021-05-18)
##  os       Ubuntu 20.04.2 LTS          
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  ctype    en_US.UTF-8                 
##  tz       UTC                         
##  date     2021-08-10                  
## 
## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
##  package     * version date       lib source        
##  bookdown      0.22    2021-04-22 [1] RSPM (R 4.1.0)
##  bslib         0.2.5.1 2021-05-18 [1] RSPM (R 4.1.0)
##  cli           3.0.1   2021-07-17 [2] RSPM (R 4.1.0)
##  digest        0.6.27  2020-10-24 [2] RSPM (R 4.1.0)
##  evaluate      0.14    2019-05-28 [2] RSPM (R 4.1.0)
##  htmltools     0.5.1.1 2021-01-22 [1] RSPM (R 4.1.0)
##  jquerylib     0.1.4   2021-04-26 [1] RSPM (R 4.1.0)
##  jsonlite      1.7.2   2020-12-09 [2] RSPM (R 4.1.0)
##  knitr         1.33    2021-04-24 [2] RSPM (R 4.1.0)
##  magrittr      2.0.1   2020-11-17 [2] RSPM (R 4.1.0)
##  R6            2.5.0   2020-10-28 [2] RSPM (R 4.1.0)
##  rlang         0.4.11  2021-04-30 [2] RSPM (R 4.1.0)
##  rmarkdown     2.10    2021-08-06 [1] RSPM (R 4.1.0)
##  sass          0.4.0   2021-05-12 [1] RSPM (R 4.1.0)
##  sessioninfo   1.1.1   2018-11-05 [2] RSPM (R 4.1.0)
##  stringi       1.7.3   2021-07-16 [2] RSPM (R 4.1.0)
##  stringr       1.4.0   2019-02-10 [2] RSPM (R 4.1.0)
##  withr         2.4.2   2021-04-18 [2] RSPM (R 4.1.0)
##  xfun          0.25    2021-08-06 [2] RSPM (R 4.1.0)
## 
## [1] /__w/_temp/Library
## [2] /usr/local/lib/R/site-library
## [3] /usr/local/lib/R/library
```

## Patrocinadores {-}

Agradecemos a nuestros patrocinadores:

<a href="https://comunidadbioinfo.github.io/es/post/cs_and_s_event_fund_award/#.YJH-wbVKj8A"><img src="https://comunidadbioinfo.github.io/post/2021-01-27-cs_and_s_event_fund_award/spanish_cs_and_s_award.jpeg" width="400px" align="center"/></a>

<a href="https://www.r-consortium.org/"><img src="https://www.r-consortium.org/wp-content/uploads/sites/13/2016/09/RConsortium_Horizontal_Pantone.png" width="400px" align="center"/></a>
