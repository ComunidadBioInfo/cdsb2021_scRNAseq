# Introducción a Seurat

Instructor: [Kevin E. Meza-Landeros](https://twitter.com/KevsGenomic)

## Diapositivas

Presentación: [aquí](https://docs.google.com/presentation/d/18ZCddwDD9lY8j4gmt1fO-Xqa8qJwCh_Zu_kXSjluW1Q/edit?usp=sharing) 

## Una perspectiva diferente 

**Seurat** es un paquete R diseñado para control de calidad, análisis y exploración de datos de secuencia de ARN de una sola célula. Seurat tiene como objetivo permitir a los usuarios identificar e interpretar fuentes de heterogeneidad a partir de mediciones transcriptómicas unicelulares e integrar diversos tipos de datos unicelulares. 

Seurat es desarrollado y mantenido por el laboratorio de [Satija](https://satijalab.org/seurat/authors.html) y se publica bajo la Licencia Pública GNU (GPL 3.0).

En este tutorial se ve como procesar los datos de scRNAseq con un nuevo paquete. Los pasos a realizar son en esencia los mismos que ya revisamos con el tutorial de la OSCA de RStudio.  
No olvides nunca que el paquete mas adecuado y que deberás utilizar dependerá mayoritariamente de tus datos y el procesamiento que se adecúe a estos.  

**Además... siempre es bueno diversos puntos de vista sobre las cosas, no es así?**

Aprende mas sobre Seurat: [aquí](https://satijalab.org/seurat/)

## Kick-start

En este tutorial partimos a partir de que ya se tienen los archivos FASTQ resultados de secuenciación.  

- ¿Con qué datos estoy trabajando?  
Peripheral Blood Mononuclear Cells **(PBMC)** disponibles gratuitamente de **10X Genomics**. Son en total 2,700 céluas únicas secuenciadas con **Illumina NextSeq 500**.
Puedes descargar los datos de [aqui](https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) (7.3MB).  
Descarga el archivo comprimido y procede a descomprimirlo. Se creara el siguiente directorio *filtered_gene_bc_matrices/hg19/*, aquí estarán los archivos que necesitaremos.

Este tutorial solo es la punta del *iceberg* de lo que se puede hacer con la paquetera de Seurat. Para comenzar a sumergirte en este mundo no dudes en visitar la página oficial mantenida por Satija Lab [Vignettes](https://satijalab.org/seurat/articles/get_started.html)

A continuación estableceremos nuestros directorio de trabajo y leeremos los datos anteriores.  
La función Read10X () lee en la salida de cellranger de 10X (de donde se obtuvieron los FASTQs), devolviendo una matriz de recuento única identificada molecularmente (UMI). Los valores en esta matriz representan el número de moléculas para cada característica (es decir, gen; fila) que se detectan en cada celda (columna).


```r
## Cargar paquetes de R
library(dplyr)
library(Seurat)
library(patchwork)
```



```r
# Load the PBMC dataset
proydir <- "/mnt/BioAdHoc/Groups/vd-vijay/kmlanderos/CDSB/clustering/"
pbmc.data <- Read10X(data.dir = paste0(proydir, "data/filtered_gene_bc_matrices/hg19/"))
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
```

Veamos la estructura del Objeto de Seurat


```r
str(pbmc)
```

- ¿Cómo se ven los datos en una matriz de recuento?  
Examinemos algunos genes en las primeras treinta células. El . los valores en la matriz representan ceros (no se detectan moléculas). Dado que la mayoría de los valores en una matriz scRNA-seq son 0, Seurat utiliza una representación de matriz dispersa (*sparse matrix*) siempre que sea posible. Esto da como resultado un ahorro significativo de memoria y velocidad.  
EN ESTE CASO UNA MATRIZ NO DISPERSA OCUPA 27 VECES MAS ESPACIO!


```r
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size
dense.size / sparse.size
```
## Quality Control

Algunas métricas de control de calidad comúnmente utilizadas por la comunidad incluyen:

- El número de genes únicos detectados en cada célula.
    - Las células de baja calidad o las gotitas vacías suelen tener muy pocos genes.
    - Los dobletes o multipletes celulares pueden exhibir un recuento de genes aberrantemente alto
- De manera similar, el número total de moléculas detectadas dentro de una célula (se correlaciona fuertemente con genes únicos)
- El porcentaje de lecturas que se asignan al genoma mitocondrial
    - Las células de baja calidad / moribundas a menudo exhiben una extensa contaminación mitocondrial
    - Calculamos métricas de control de calidad mitocondrial con la función PercentageFeatureSet (), que calcula el porcentaje de recuentos que se originan a partir de un conjunto de características.

El operador [[puede agregar columnas a los metadatos del objeto. Este es un gran lugar para almacenar estadísticas de control de calidad. Entonces calculamos y añadimos la cantidad de lecturas que corresponden al genoma mitocondrial.


```r
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
```
Visualizamos las métricas de control de calidad mencionadas anteriormente como un diagrama de violín. Ademas vemos la correlacion entre el numero de moleculaas de RNA detectadas en cada clula con el número de genes únicos y con el porcentaje de lecturas que corresponden a mtADN.


```r
plot <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
print(plot)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)
```
Finalmente filtramos aquellas células que se salen de los estndares de cada uno de los parámetros.


```r
# Filter
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
```
- Dónde se almacenan la métricas de QC en Seurat?
Estan almacenadas en la seccion de **meta-data** del objeto Seurat


```r
head(pbmc@meta.data, 5)
```

## Normalización

De forma predeterminada, se emplea un método de normalización de escala global **"LogNormalize"** que normaliza las medidas de expresión de características para cada celda por la expresión total, multiplica esto por un factor de escala (10.000 por defecto) y transforma el resultado en logaritmos. Los valores normalizados se almacenan en pbmc [["RNA"]] @ data 


```r
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
```

## Detección de genes (caractersticas) altamente variables

A continuación, calculamos un subconjunto de **características que exhiben una alta variación de célula a célula en el conjunto de datos** (es decir, están altamente expresadas en algunas células y poco expresadas en otras). El equipo de Seurat y otros equipos han descubierto que centrarse en estos genes en el análisis posterior ayuda a resaltar la señal biológica en conjuntos de datos unicelulares.

Nuestro procedimiento en Seurat se describe en detalle aquí y mejora las versiones anteriores al modelar directamente la relación de varianza media inherente a los datos de una sola celda, y se implementa en la función *FindVariableFeatures()*. De forma predeterminada, **devolvemos 2000 características por conjunto de datos** (aunque se puede modificar). Estos se utilizarán en análisis posteriores, como PCA. 


```r
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
print(plot1 + plot2)
```
## Escalar los datos

A continuación, aplicamos una transformación lineal ("escalado") que es un paso de preprocesamiento estándar antes de las técnicas de reducción dimensional como PCA. La función *ScaleData()*:
- Cambia la expresión de cada gen, de modo que la expresión media en las células sea 0
- Escala la expresión de cada gen, de modo que la varianza entre las células sea 1
     - Este paso otorga el mismo peso en los análisis posteriores, de modo que los genes altamente expresados no dominen
Los resultados de esto se almacenan en pbmc [["RNA"]] @ scale.data


```r
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
```

## Reducción dimensional lineal

A continuación, realizamos PCA sobre los datos escalados. De forma predeterminada, solo las características variables determinadas previamente se utilizan como entrada, pero se pueden definir mediante el argumento de características si desea elegir un subconjunto diferente. 


```r
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
```

Seurat proporciona varias formas útiles de visualizar tanto las celdas como las características que definen el PCA, incluidas *VizDimReduction()*, *DimPlot()* y *DimHeatmap()* 

Puedes examinar y visualice los resultados de PCA de diferentes formas 


```r
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

p1 <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
p2 <- DimPlot(pbmc, reduction = "pca")
print(p1)
print(p2)
```

En particular, *DimHeatmap()* permite una fácil exploración de las fuentes primarias de heterogeneidad en un conjunto de datos y puede ser útil cuando se intenta decidir qué PC incluir para análisis posteriores posteriores. Tanto las celdas como las características se ordenan de acuerdo con sus puntajes de PCA. Establecer celdas en un número traza las celdas "extremas" en ambos extremos del espectro, lo que acelera drásticamente el trazado de grandes conjuntos de datos. Aunque claramente es un análisis supervisado, consideramos que esta es una herramienta valiosa para explorar conjuntos de características correlacionadas. 


```r
p3 <- DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
p4 <- DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
print(p3)
print(p4)
```

## Determinar la dimensionalidad del conjunto de datos 

Para superar el extenso ruido técnico en cualquier característica única para los datos de scRNA-seq, Seurat agrupa las células en función de sus puntuaciones de PCA, y cada PC representa esencialmente una "metafunción" que combina información en un conjunto de características correlacionadas. Por lo tanto, los componentes principales principales representan una compresión sólida del conjunto de datos. **Sin embargo, ¿cuántos componentes deberíamos elegir incluir? 10? 20? 100?**

En Macosko et al, implementamos una prueba de remuestreo inspirada en el **procedimiento JackStraw**. Permutamos aleatoriamente un subconjunto de los datos (1% por defecto) y volvemos a ejecutar PCA, construyendo una "distribución nula" de puntuaciones de características, y repetimos este procedimiento. Identificamos PC "importantes" como aquellas que tienen un gran enriquecimiento de características de bajo valor p. 


```r
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
```

La función **JackStrawPlot()** proporciona una herramienta de visualización para comparar la distribución de los valores p para cada PC con una distribución uniforme (línea discontinua). Las PC "significativas" mostrarán un gran enriquecimiento de funciones con valores p bajos (curva sólida por encima de la línea discontinua). En este caso, parece que hay una fuerte caída en la importancia después de los primeros 10-12 PCs. 


```r
p1 <- JackStrawPlot(pbmc, dims = 1:15)
```

Un método heurístico alternativo genera un **"diagrama de codo (Elbow Plot)"**: una clasificación de componentes principales basada en el porcentaje de varianza explicada por cada uno (función ElbowPlot ()). En este ejemplo, podemos observar un "codo" alrededor de PC9-10, lo que sugiere que la mayor parte de la señal verdadera se captura en las primeras 10 PC. 


```r
p2 <- ElbowPlot(pbmc)
print(p1)
print(p2)
```

## Clustering 


```r
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
```

## Reducción dimensional no lineal (UMAP/tSNE)

Seurat ofrece varias técnicas de reducción dimensional no lineal, como **tSNE** y **UMAP**, para visualizar y explorar estos conjuntos de datos. El objetivo de estos algoritmos es aprender la variedad subyacente de los datos para colocar celdas similares juntas en un espacio de baja dimensión. Las celdas dentro de los grupos basados en gráficos determinados anteriormente deben ubicarse conjuntamente en estos gráficos de reducción de dimensión. Como entrada para UMAP y tSNE, sugerimos usar las mismas PC como entrada para el análisis de agrupamiento. 


```r
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
plot <- DimPlot(pbmc, reduction = "umap")
print(plot)
```

Puede guardar el objeto en este punto para que se pueda volver a cargar fácilmente sin tener que volver a ejecutar los pasos computacionalmente intensivos realizados anteriormente o compartir fácilmente con los colaboradores. 


```r
saveRDS(pbmc, file = "output/pbmc_tutorial.rds")
```

## Caracteristicas diferencialmente expresadas (biomarcadores de los clusters)

Seurat puede ayudarlo a encontrar marcadores que definan clústeres mediante expresión diferencial. De forma predeterminada, identifica **marcadores positivos y negativos de un solo grupo (especificado en ident.1), en comparación con todas las demás células**. FindAllMarkers () automatiza este proceso para todos los clústeres, pero **también se pueden comparar grupos de clústeres entre sí o contra todas las celdas**.

El argumento min.pct requiere que se detecte una característica en un porcentaje mínimo en cualquiera de los dos grupos de celdas, y el argumento thresh.test requiere que una característica se exprese diferencialmente (en promedio) en alguna cantidad entre los dos grupos. Puede establecer ambos en 0, pero con un aumento dramático en el tiempo, ya que esto probará una gran cantidad de características que probablemente no sean altamente discriminatorias. Como otra opción para acelerar estos cálculos, se puede configurar el número máximo de celdas por identificador. Esto reducirá la resolución de cada clase de identidad para que no tenga más celdas que las que se establezcan. Si bien generalmente habrá una pérdida de potencia, los aumentos de velocidad pueden ser significativos y es probable que las características expresadas de manera más diferencial aún se eleven a la cima. 


```r
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
```

Seurat tiene **varias pruebas de expresión diferencial** que se pueden configurar con el parámetro test.use (consulte nuestra viñeta DE para obtener más detalles). Por ejemplo, la **prueba ROC** devuelve el "poder de clasificación" para cualquier marcador individual (que varía de 0 - aleatorio a 1 - perfecto) .


```r
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
```

Se incluyen varias herramientas para visualizar la expresión de los marcadores. VlnPlot () (muestra distribuciones de probabilidad de expresión entre clústeres) y FeaturePlot () (visualiza la expresión de características en un gráfico tSNE o PCA) son nuestras visualizaciones más utilizadas. También sugerimos explorar RidgePlot (), CellScatter () y DotPlot () como métodos adicionales para ver su conjunto de datos. 


```r
p1 <- VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

## you can plot raw counts as well
p2 <- VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

p3 <- FeaturePlot(pbmc, features = c(
    "MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"
))

# DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
p4 <- DoHeatmap(pbmc, features = top10$gene) + NoLegend()
print(p1)
print(p2)
print(p3)
print(p4)
```

## Assigning cell type identity to clusters

Podemos usar marcadores canónicos para hacer coincidir fácilmente la agrupación imparcial con los tipos de células conocidos.


```r
new.cluster.ids <- c(
    "Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet"
)
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

p1 <- DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
print(p1)
```

## Guardar Resultados


```r
saveRDS(pbmc, file = "output/pbmc3k_final.rds")
```

## Detalles de la sesión de R


```r
## Información de la sesión de R
Sys.time()
```

```
## [1] "2021-08-09 03:55:20 UTC"
```

```r
proc.time()
```

```
##    user  system elapsed 
##   3.711   0.235   3.857
```

```r
options(width = 120)
sessioninfo::session_info()
```

```
## Registered S3 method overwritten by 'cli':
##   method     from         
##   print.boxx spatstat.geom
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
##  date     2021-08-09                  
## 
## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
##  package         * version date       lib source        
##  abind             1.4-5   2016-07-21 [1] RSPM (R 4.1.0)
##  assertthat        0.2.1   2019-03-21 [1] RSPM (R 4.1.0)
##  bookdown          0.22    2021-04-22 [1] RSPM (R 4.1.0)
##  bslib             0.2.5.1 2021-05-18 [1] RSPM (R 4.1.0)
##  cli               3.0.1   2021-07-17 [2] RSPM (R 4.1.0)
##  cluster           2.1.2   2021-04-17 [3] CRAN (R 4.1.0)
##  codetools         0.2-18  2020-11-04 [3] CRAN (R 4.1.0)
##  colorspace        2.0-2   2021-06-24 [1] RSPM (R 4.1.0)
##  cowplot           1.1.1   2020-12-30 [1] RSPM (R 4.1.0)
##  crayon            1.4.1   2021-02-08 [2] RSPM (R 4.1.0)
##  data.table        1.14.0  2021-02-21 [1] RSPM (R 4.1.0)
##  DBI               1.1.1   2021-01-15 [1] RSPM (R 4.1.0)
##  deldir            0.2-10  2021-02-16 [1] RSPM (R 4.1.0)
##  digest            0.6.27  2020-10-24 [2] RSPM (R 4.1.0)
##  dplyr           * 1.0.7   2021-06-18 [1] RSPM (R 4.1.0)
##  ellipsis          0.3.2   2021-04-29 [2] RSPM (R 4.1.0)
##  evaluate          0.14    2019-05-28 [2] RSPM (R 4.1.0)
##  fansi             0.5.0   2021-05-25 [2] RSPM (R 4.1.0)
##  fastmap           1.1.0   2021-01-25 [2] RSPM (R 4.1.0)
##  fitdistrplus      1.1-5   2021-05-28 [1] RSPM (R 4.1.0)
##  future            1.21.0  2020-12-10 [1] RSPM (R 4.1.0)
##  future.apply      1.7.0   2021-01-04 [1] RSPM (R 4.1.0)
##  generics          0.1.0   2020-10-31 [1] RSPM (R 4.1.0)
##  ggplot2           3.3.5   2021-06-25 [1] RSPM (R 4.1.0)
##  ggrepel           0.9.1   2021-01-15 [1] RSPM (R 4.1.0)
##  ggridges          0.5.3   2021-01-08 [1] RSPM (R 4.1.0)
##  globals           0.14.0  2020-11-22 [1] RSPM (R 4.1.0)
##  glue              1.4.2   2020-08-27 [2] RSPM (R 4.1.0)
##  goftest           1.2-2   2019-12-02 [1] RSPM (R 4.1.0)
##  gridExtra         2.3     2017-09-09 [1] RSPM (R 4.1.0)
##  gtable            0.3.0   2019-03-25 [1] RSPM (R 4.1.0)
##  htmltools         0.5.1.1 2021-01-22 [1] RSPM (R 4.1.0)
##  htmlwidgets       1.5.3   2020-12-10 [1] RSPM (R 4.1.0)
##  httpuv            1.6.1   2021-05-07 [1] RSPM (R 4.1.0)
##  httr              1.4.2   2020-07-20 [2] RSPM (R 4.1.0)
##  ica               1.0-2   2018-05-24 [1] RSPM (R 4.1.0)
##  igraph            1.2.6   2020-10-06 [1] RSPM (R 4.1.0)
##  irlba             2.3.3   2019-02-05 [1] RSPM (R 4.1.0)
##  jquerylib         0.1.4   2021-04-26 [1] RSPM (R 4.1.0)
##  jsonlite          1.7.2   2020-12-09 [2] RSPM (R 4.1.0)
##  KernSmooth        2.23-20 2021-05-03 [3] CRAN (R 4.1.0)
##  knitr             1.33    2021-04-24 [2] RSPM (R 4.1.0)
##  later             1.2.0   2021-04-23 [1] RSPM (R 4.1.0)
##  lattice           0.20-44 2021-05-02 [3] CRAN (R 4.1.0)
##  lazyeval          0.2.2   2019-03-15 [1] RSPM (R 4.1.0)
##  leiden            0.3.9   2021-07-27 [1] RSPM (R 4.1.0)
##  lifecycle         1.0.0   2021-02-15 [2] RSPM (R 4.1.0)
##  listenv           0.8.0   2019-12-05 [1] RSPM (R 4.1.0)
##  lmtest            0.9-38  2020-09-09 [1] RSPM (R 4.1.0)
##  magrittr          2.0.1   2020-11-17 [2] RSPM (R 4.1.0)
##  MASS              7.3-54  2021-05-03 [3] CRAN (R 4.1.0)
##  Matrix            1.3-4   2021-06-01 [3] RSPM (R 4.1.0)
##  matrixStats       0.60.0  2021-07-26 [1] RSPM (R 4.1.0)
##  mgcv              1.8-36  2021-06-01 [3] RSPM (R 4.1.0)
##  mime              0.11    2021-06-23 [2] RSPM (R 4.1.0)
##  miniUI            0.1.1.1 2018-05-18 [1] RSPM (R 4.1.0)
##  munsell           0.5.0   2018-06-12 [1] RSPM (R 4.1.0)
##  nlme              3.1-152 2021-02-04 [3] CRAN (R 4.1.0)
##  parallelly        1.27.0  2021-07-19 [1] RSPM (R 4.1.0)
##  patchwork       * 1.1.1   2020-12-17 [1] RSPM (R 4.1.0)
##  pbapply           1.4-3   2020-08-18 [1] RSPM (R 4.1.0)
##  pillar            1.6.2   2021-07-29 [2] RSPM (R 4.1.0)
##  pkgconfig         2.0.3   2019-09-22 [2] RSPM (R 4.1.0)
##  plotly            4.9.4.1 2021-06-18 [1] RSPM (R 4.1.0)
##  plyr              1.8.6   2020-03-03 [1] RSPM (R 4.1.0)
##  png               0.1-7   2013-12-03 [1] RSPM (R 4.1.0)
##  polyclip          1.10-0  2019-03-14 [1] RSPM (R 4.1.0)
##  promises          1.2.0.1 2021-02-11 [1] RSPM (R 4.1.0)
##  purrr             0.3.4   2020-04-17 [2] RSPM (R 4.1.0)
##  R6                2.5.0   2020-10-28 [2] RSPM (R 4.1.0)
##  RANN              2.6.1   2019-01-08 [1] RSPM (R 4.1.0)
##  RColorBrewer      1.1-2   2014-12-07 [1] RSPM (R 4.1.0)
##  Rcpp              1.0.7   2021-07-07 [2] RSPM (R 4.1.0)
##  RcppAnnoy         0.0.19  2021-07-30 [1] RSPM (R 4.1.0)
##  reshape2          1.4.4   2020-04-09 [1] RSPM (R 4.1.0)
##  reticulate        1.20    2021-05-03 [1] RSPM (R 4.1.0)
##  rlang             0.4.11  2021-04-30 [2] RSPM (R 4.1.0)
##  rmarkdown         2.9     2021-06-15 [1] RSPM (R 4.1.0)
##  ROCR              1.0-11  2020-05-02 [1] RSPM (R 4.1.0)
##  rpart             4.1-15  2019-04-12 [3] CRAN (R 4.1.0)
##  Rtsne             0.15    2018-11-10 [1] RSPM (R 4.1.0)
##  sass              0.4.0   2021-05-12 [1] RSPM (R 4.1.0)
##  scales            1.1.1   2020-05-11 [1] RSPM (R 4.1.0)
##  scattermore       0.7     2020-11-24 [1] RSPM (R 4.1.0)
##  sctransform       0.3.2   2020-12-16 [1] RSPM (R 4.1.0)
##  sessioninfo       1.1.1   2018-11-05 [2] RSPM (R 4.1.0)
##  Seurat          * 4.0.3   2021-06-10 [1] RSPM (R 4.1.0)
##  SeuratObject    * 4.0.2   2021-06-09 [1] RSPM (R 4.1.0)
##  shiny             1.6.0   2021-01-25 [1] RSPM (R 4.1.0)
##  spatstat.core     2.3-0   2021-07-16 [1] RSPM (R 4.1.0)
##  spatstat.data     2.1-0   2021-03-21 [1] RSPM (R 4.1.0)
##  spatstat.geom     2.2-2   2021-07-12 [1] RSPM (R 4.1.0)
##  spatstat.sparse   2.0-0   2021-03-16 [1] RSPM (R 4.1.0)
##  spatstat.utils    2.2-0   2021-06-14 [1] RSPM (R 4.1.0)
##  stringi           1.7.3   2021-07-16 [2] RSPM (R 4.1.0)
##  stringr           1.4.0   2019-02-10 [2] RSPM (R 4.1.0)
##  survival          3.2-11  2021-04-26 [3] CRAN (R 4.1.0)
##  tensor            1.5     2012-05-05 [1] RSPM (R 4.1.0)
##  tibble            3.1.3   2021-07-23 [2] RSPM (R 4.1.0)
##  tidyr             1.1.3   2021-03-03 [1] RSPM (R 4.1.0)
##  tidyselect        1.1.1   2021-04-30 [1] RSPM (R 4.1.0)
##  utf8              1.2.2   2021-07-24 [2] RSPM (R 4.1.0)
##  uwot              0.1.10  2020-12-15 [1] RSPM (R 4.1.0)
##  vctrs             0.3.8   2021-04-29 [2] RSPM (R 4.1.0)
##  viridisLite       0.4.0   2021-04-13 [1] RSPM (R 4.1.0)
##  withr             2.4.2   2021-04-18 [2] RSPM (R 4.1.0)
##  xfun              0.24    2021-06-15 [2] RSPM (R 4.1.0)
##  xtable            1.8-4   2019-04-21 [1] RSPM (R 4.1.0)
##  zoo               1.8-9   2021-03-09 [1] RSPM (R 4.1.0)
## 
## [1] /__w/_temp/Library
## [2] /usr/local/lib/R/site-library
## [3] /usr/local/lib/R/library
```

## Patrocinadores {-}

Agradecemos a nuestros patrocinadores:

<a href="https://comunidadbioinfo.github.io/es/post/cs_and_s_event_fund_award/#.YJH-wbVKj8A"><img src="https://comunidadbioinfo.github.io/post/2021-01-27-cs_and_s_event_fund_award/spanish_cs_and_s_award.jpeg" width="400px" align="center"/></a>

<a href="https://www.r-consortium.org/"><img src="https://www.r-consortium.org/wp-content/uploads/sites/13/2016/09/RConsortium_Horizontal_Pantone.png" width="400px" align="center"/></a>
