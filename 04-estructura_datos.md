# Estructura e importe de datos

Instructoras: [Elisa Márquez Zavala](https://twitter.com/naielisha), [Citlali Gil Aguillon](http://twitter.com/argininaa)

Contenido adaptado de [CDSB2020: Introducción a scRNA-seq, estructura e importe de datos
](https://comunidadbioinfo.github.io/cdsb2020/scRNAseq/02-data-infrastructure-and-import.html#1) de [Leonardo Collado Torres](https://github.com/lcolladotor) 

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
## [1] "2021-08-10 07:24:54 UTC"
```

```r
proc.time()
```

```
##    user  system elapsed 
##   0.413   0.130   0.437
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
