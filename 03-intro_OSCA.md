# Introducción a RNA-seq de célula única (scRNA-seq) con Bioconductor y al libro de OSCA

Instructoras: [Elisa Márquez Zavala](https://twitter.com/naielisha), [Citlali Gil Aguillon](http://twitter.com/argininaa)

Contenido adaptado del [Curso de RNASeq](https://lcolladotor.github.io/rnaseq_LCG-UNAM_2021/index.html) de [Leonardo Collado Torres](https://github.com/lcolladotor) y de 

<iframe width="560" height="315" src="https://www.youtube.com/embed/c0Ch_sXiGDQ" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

[_Original Notes in English_](https://docs.google.com/document/d/1VZ4dOesHjbWhrvSIXO8Ne1Qk4ZdGEhAySNviZOHPqjI/edit?usp=sharing)

## Bioconductor


* CRAN, the R package repository: https://cran.r-project.org/
* CRAN task views: https://cran.r-project.org/web/views/


"Bioconductor proporciona herramientas para el análisis y la comprensión de datos genómicos de alto rendimiento. Bioconductor utiliza el lenguaje de programación estadístico R y es de código abierto y desarrollo abierto. Tiene dos lanzamientos cada año y una comunidad de usuarios activa. Bioconductor también está disponible como AMI (Imagen de máquina de Amazon) e imágenes de Docker." 

* https://www.bioconductor.org/

* Where do I start using Bioconductor? http://lcolladotor.github.io/2014/10/16/startbioc/#.XqxNGRNKiuo

Básicamente es un repositorio con reglas o estándares  para el análisis y la comprensión de datos genómicos de alto rendimiento.

Para conocer sobre Bioconductor podemos ir a:
https://www.bioconductor.org/ y dar click en *About*


### Equipos y consejos 

Es conformado por diversos equipos y consejos ([Asesores científicos](https://www.bioconductor.org/about/scientific-advisory-board/), [técnicos](https://www.bioconductor.org/about/technical-advisory-board/) y de la [comunidad](https://www.bioconductor.org/about/community-advisory-board/)). Por ejemplo Leonardo Collado.

  - Científicos : Proporciona orientación externa y supervisión de la dirección científica del proyecto y está compuesto por líderes en el análisis estadístico de datos genómicos. 
  - Técnicos: Desarrollar estrategias para asegurar que la parte técnica de la infraestructura sea apropiada a largo plazo (manejo de paquetes, sitio web, slack, etc)
  - Comunidad : Empoderar a las comunidades de usuarios y desarrolladores mediante la coordinación de actividades de capacitación y divulgación.

Dentro del equipo *core* que mantiene a Bioconductor y apoya con las dudas (https://www.bioconductor.org/about/core-team/) hay gente a la que Bioconductor le paga por mantener los repositorios, lo cual lo hace diferente de CRAN. El tener gente que oficialmente sabe cómo ayudarte y tiene el tiempo para hacerlo crea una mejor experiencia para los usuarios y los desarrolladores.


### Encontrando paquetes de Bioconductor 

* Tipos de paquetes
Hay 4 tipos de paquetes que aceptan

  * [Software](http://bioconductor.org/packages/release/bioc/): tipo principal de paquete BioC, en su mayoría aportado por el usuario. Es un paquete con un tipo de análisis específico. Algunos los hacen gente pagada directamente por Bioconductor

  * [Annotation]( http://bioconductor.org/packages/release/data/annotation/): facilita la interacción con bases de datos genómicas muy utilizadas

  * [Experiment](http://bioconductor.org/packages/release/data/experiment/): contienen datos para algún artículo o datos que se usan en ejemplos más exhaustivos, en su mayoría aportados por el usuario. ~<5 Mb

  * [Workflows](http://bioconductor.org/packages/release/workflows/): demuestran como puedes usar varios paquetes de Bioconductor para ciertos tipos de análisis

* Para descubrir paquetes:
  - Software: http://bioconductor.org/packages/release/bioc/
  - Annotation: http://bioconductor.org/packages/release/data/annotation/
  - Experiment Data: http://bioconductor.org/packages/release/data/experiment/
  - Workflows: http://bioconductor.org/packages/release/workflows/


Las listas de cada tipo de paquete se ven algo así:

| Package            | Maintainer              | Title                                               |
|--------------------|-------------------------|-----------------------------------------------------|
| Nombre del paquete | Quién lo mantiene       | Título completo                                     |
| recount3           | Leonardo Collado-Torres | Explore and download data from the recount3 project |


* Paquetes de R de Leo: https://lcolladotor.github.io/pkgs/

Sin embargo, estas listas no son muy amigables si queremos explorar por lo que podemos usar `biocViews`


* Encontrando paquetes a través de `biocViews`: http://bioconductor.org/packages/release/BiocViews.html#___Software
    - Estructura tipo árbol
    - Son 4 árboles principales: software, annotation, experiment, workflow
    - Dentro de cada árbol, un paquete puede ser parte de varias ramas, por ejemplo, [recount3](https://www.bioconductor.org/packages/release/bioc/html/recount3.html) está dentro de todas estas ramas:
      - <span style="color:Purple">Software</span>
        - <span style="color:Purple">AssayDomain</span>
          - <span style="color:Purple">GeneExpression</span>
        - <span style="color:Purple">BiologicalQuestion</span>
          - <span style="color:Purple">DifferentialExpression</span>
          - <span style="color:Purple">Coverage</span>
        - <span style="color:Purple">Infrastructure</span>
          - <span style="color:Purple">DataImport</span>
        - <span style="color:Purple">Technology</span>
          - <span style="color:Purple">Sequencing</span>
            - <span style="color:Purple">RNASeq</span>
    - Tiene una búsqueda de texto simple
    - Ejemplo: Software → WorkflowStep → Visualization → http://bioconductor.org/packages/release/BiocViews.html#___Visualization (486 paquetes en BioC 3.11 abril-octubre 2020, 506 en BioC 3.12 octubre 2020-abril 2021, 529 en BioC 3.13 agosto 2021)

### Estructura de un paquete de BioC

* Usa `https://bioconductor.org/packages/<pkg_name>`
    - Ejemplo: https://bioconductor.org/packages/recount
    - Otro ejemplo: https://bioconductor.org/packages/SummarizedExperiment
* Badges (etiquetas): rápidamente podemos evaluar como está

![](https://www.bioconductor.org/shields/availability/release/SummarizedExperiment.svg) : ¿En qué plataformas funciona?

![](https://www.bioconductor.org/shields/downloads/release/SummarizedExperiment.svg) : ¿Qué tan descargado es?

![](https://www.bioconductor.org/shields/posts/SummarizedExperiment.svg) : ¿Se han hecho preguntas del paquete en los últimos 6 meses? (respondidas/hechas)

![](https://www.bioconductor.org/shields/years-in-bioc/SummarizedExperiment.svg) : ¿Cuánto tiempo lleva en Bioconductor?

![](https://www.bioconductor.org/shields/build/release/bioc/SummarizedExperiment.svg) : ¿Funciona en las máquinas de bioconductor?

![](https://www.bioconductor.org/shields/lastcommit/release/bioc/SummarizedExperiment.svg) : ¿Cuándo fue la última vez que lo actualizaron?

![](https://www.bioconductor.org//shields/dependencies/release/SummarizedExperiment.svg) : Número de dependencias recursivas necesarias para instalar el paquete 

* Parráfo de descripción del paquete
* Cómo citar al paquete de Bioconductor
* Cómo instalarlo. Más detalles en http://bioconductor.org/install/
* Documentación
    - Una líga por cada vignette en formato PDF o HTML. Es la documentación **principal**!
    - Una vignette es donde lxs autores del paquete explican cómo usar las diferentes funciones del paquete y en qué orden
* Detalles
    - Términos de `biocViews`
    - Cómo se relaciona a otros paquetes (depends, imports, linking to, suggests, depends on me, ...)
    - URL: donde puedes encontrar el código fuente (nos puede dar más infor)
    - BugReports: donde puedes pedir ayuda
* Más detalles sobre el paquete
    - Estadísticas de descargas

### Las dos ramas de Bioconductor: release y devel

* Dos ramas
    - `release`, actualmente 3.13
    - `devel`, actualmente 3.14
    
    - Bioconductor version 3.14 (Development)
 https://bioconductor.org/packages/devel/BiocViews.html#___Software
    - Ejemplo: http://bioconductor.org/packages/devel/bioc/html/recount.html
* Bioconductor tiene es actualizado cada 6 meses (abril y octubre). R lo actualizan 1 vez al año (abril).
* Todo el software lo prueban en macOS, Windows y linux
    - Ejemplo: http://bioconductor.org/checkResults/release/bioc-LATEST/recount/ y http://bioconductor.org/checkResults/devel/bioc-LATEST/recount/
* Resumen BioC 3.13 http://bioconductor.org/news/bioc_3_13_release/
    - Blog post en LIBD rstats club: Quick overview on the new Bioconductor 3.8 release http://research.libd.org/rstatsclub/2018/11/02/quick-overview-on-the-new-bioconductor-3-8-release/


### Cursos y eventos

* http://bioconductor.org/help/events/
* http://bioconductor.org/help/course-materials/
* BioC2021: conferencia principal anual https://bioc2021.bioconductor.org/
    - Talleres del BioC2019: https://rebrand.ly/biocworkshops2019

<blockquote class="twitter-tweet"><p lang="en" dir="ltr">Teach online data science, bioinformatics, or other computational skills interactively using the Orchestra platform:<a href="https://t.co/r4aJ2xAZbh">https://t.co/r4aJ2xAZbh</a> <br>Nearly 50 workshop environments preloaded with <a href="https://twitter.com/hashtag/jupyter?src=hash&amp;ref_src=twsrc%5Etfw">#jupyter</a>, <a href="https://twitter.com/hashtag/rstudio?src=hash&amp;ref_src=twsrc%5Etfw">#rstudio</a>, <a href="https://twitter.com/hashtag/shell?src=hash&amp;ref_src=twsrc%5Etfw">#shell</a>. <a href="https://twitter.com/hashtag/rstats?src=hash&amp;ref_src=twsrc%5Etfw">#rstats</a>, or <a href="https://twitter.com/hashtag/python?src=hash&amp;ref_src=twsrc%5Etfw">#python</a>.<a href="https://twitter.com/NIHSTRIDES?ref_src=twsrc%5Etfw">@NIHSTRIDES</a> <a href="https://twitter.com/NIHDataScience?ref_src=twsrc%5Etfw">@NIHDataScience</a> <a href="https://twitter.com/Bioconductor?ref_src=twsrc%5Etfw">@Bioconductor</a> <a href="https://t.co/HyWVLBJxGU">pic.twitter.com/HyWVLBJxGU</a></p>&mdash; Sean Davis (@seandavis12) <a href="https://twitter.com/seandavis12/status/1348291911090626560?ref_src=twsrc%5Etfw">January 10, 2021</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

* Talleres de la CDSB, como los talleres de CDSB 2020: https://comunidadbioinfo.github.io/post/cdsb2020-building-workflows-with-rstudio-and-scrnaseq-with-bioconductor/#.XmJT-Z-YU1I



### Comunidad

* Slack: https://bioc-community.herokuapp.com/
* Sitio web de ayuda: https://support.bioconductor.org/
    - Usa la(s) etiqueta(s) adecuada(s) para que lxs autores de los paquetes reciban email de forma automática
    - Pueden revisar ese sitio web y usarlo para aprender cómo en https://lcolladotor.github.io/bioc_team_ds/helping-others.html#bioconductor-support-practice-grounds
* Twitter: https://twitter.com/bioconductor
  
## Introducción a RNA-seq de célula única (scRNA-seq) con Bioconductor y al libro de OSCA  

[link a diapositivas](https://docs.google.com/presentation/d/1Jr_CtxAX4MviZF6Hb4taPOPWQl1HZW942M_pH51B4Vg/edit?usp=sharing)


## Detalles de la sesión de R


```r
## Información de la sesión de R
Sys.time()
```

```
## [1] "2021-08-12 07:39:05 UTC"
```

```r
proc.time()
```

```
##    user  system elapsed 
##   0.409   0.120   0.408
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
##  date     2021-08-12                  
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
