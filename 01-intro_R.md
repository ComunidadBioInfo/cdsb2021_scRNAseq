# Inauguración CDSB2021 {-}

<script async class="speakerdeck-embed" data-id="2183a59718804e359fedef8de75606d2" data-ratio="1.77777777777778" src="//speakerdeck.com/assets/embed.js"></script>

# Introducción a R y RStudio

Instructor: [Leonardo Collado Torres](https://comunidadbioinfo.github.io/es/authors/lcollado/)

## R

* R: es gratis, de acceso libre, utilizado para muchos campos de trabajo, fuerte en la bioinformática a través de Bioconductor
* Instalación a través de CRAN: https://cran.r-project.org/
* Para explorar que se puede hacer con R:
  - R Weekly https://rweekly.org/
  - R Bloggers https://www.r-bloggers.com/
  - Twitter https://twitter.com/search?q=%23rstats&src=typed_query
  - Twitter en español https://twitter.com/search?q=%23rstatsES&src=typed_query
  - TidyTuesday https://twitter.com/search?q=%23TidyTuesday&src=typed_query
  - DatosDeMiercoles https://twitter.com/search?q=%23datosdemiercoles&src=typed_query
* Para pedir ayuda hay muchas opciones
  - https://lcolladotor.github.io/bioc_team_ds/how-to-ask-for-help.html
  

* Material en el que estoy involucrado:
  - https://twitter.com/lcolladotor
  - https://www.youtube.com/c/LeonardoColladoTorres/playlists
  - LIBD rstats club https://docs.google.com/spreadsheets/d/1is8dZSd0FZ9Qi1Zvq1uRhm-P1McnJRd_zxdAfCRoMfA/edit?usp=sharing
  - https://twitter.com/CDSBMexico, https://twitter.com/LIBDrstats, https://twitter.com/Bioconductor
  - https://comunidadbioinfo.github.io/
  - YouTube CDSB: https://www.youtube.com/channel/UCHCdYfAXVzJIUkMoMSGiZMw
  - _Orchestrating Spatially Resolved Transcriptomics Analysis with Bioconductor_ https://lmweber.org/OSTA-book/

<img src="https://raw.githubusercontent.com/lcolladotor/rnaseq_LCG-UNAM_2021/master/icon_192.png">

<blockquote class="twitter-tweet"><p lang="es" dir="ltr">Thank you <a href="https://twitter.com/Bioconductor?ref_src=twsrc%5Etfw">@Bioconductor</a> for welcoming me &amp; providing me a foundation &amp; platform for my career<br><br>Gracias BioC por darme una oportunidad y plataforma para desarrollar mi carrera<br><br>Time to pass it on/Toca enseñar y crecer la bola de nieve ❄️ <a href="https://twitter.com/CDSBMexico?ref_src=twsrc%5Etfw">@CDSBMexico</a> 🇲🇽<a href="https://twitter.com/hashtag/BioC2021?src=hash&amp;ref_src=twsrc%5Etfw">#BioC2021</a> <a href="https://twitter.com/hashtag/rstats?src=hash&amp;ref_src=twsrc%5Etfw">#rstats</a> <a href="https://twitter.com/hashtag/rstatsES?src=hash&amp;ref_src=twsrc%5Etfw">#rstatsES</a> <a href="https://t.co/8cDXP4Pf2W">https://t.co/8cDXP4Pf2W</a> <a href="https://t.co/rpQgH8UsWW">pic.twitter.com/rpQgH8UsWW</a></p>&mdash; 🇲🇽 Leonardo Collado-Torres (@lcolladotor) <a href="https://twitter.com/lcolladotor/status/1423385604314865664?ref_src=twsrc%5Etfw">August 5, 2021</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script>

## GitHub

* Permite compartir código
* Se complementa con Git que es para tener un control de versiones de tu código
  - https://github.com/ComunidadBioInfo/cdsb2020/blob/master/presentaciones_flujos-de-trabajo/Introduccion-al-flujo-de-trabajo-orientado-a-proyectos.pdf
* Puedes tener páginas web estáticas
  - https://pages.github.com/
  - https://github.com/ComunidadBioInfo/cdsb2021_scRNAseq/. En especial https://github.com/ComunidadBioInfo/cdsb2021_scRNAseq/tree/gh-pages se convierte en https://comunidadbioinfo.github.io/cdsb2021_scRNAseq/
  - Página personal: https://github.com/lcolladotor/lcolladotor.github.com se convierte en http://lcolladotor.github.io/. Está todo hecho con https://github.com/lcolladotor/lcolladotorsource
* Tip: usen el mismo nombre de usuario en GitHub, Twitter, Gmail, etc. 
  - How to be a Modern Scientist: https://lcolladotor.github.io/bioc_team_ds/how-to-be-a-modern-scientist.html
  
## RStudio

* RStudio Desktop es gratis http://www.rstudio.com/products/rstudio/download/preview/
* Nos ayuda a realizar muchas cosas con R de forma más rápida
  - Demo `rsthemes`
  

```r
remotes::install_github(c(
    "gadenbuie/rsthemes"
))
remotes::install_cran("suncalc")
rsthemes::install_rsthemes(include_base16 = TRUE)
```

  

```r
usethis::edit_r_profile()

## From https://www.garrickadenbuie.com/project/rsthemes/
if (interactive() && requireNamespace("rsthemes", quietly = TRUE)) {
    # Set preferred themes if not handled elsewhere..
    rsthemes::set_theme_light("Solarized Light {rsthemes}") # light theme
    rsthemes::set_theme_dark("base16 Monokai {rsthemes}") # dark theme
    rsthemes::set_theme_favorite(c(
        "Solarized Light {rsthemes}",
        "base16 Monokai {rsthemes}",
        "One Dark {rsthemes}"
    ))
    # Whenever the R session restarts inside RStudio...
    setHook("rstudio.sessionInit", function(isNewSession) {
        # Automatically choose the correct theme based on time of day
        ## Used rsthemes::geolocate() once
        rsthemes::use_theme_auto(lat = 39.2891, lon = -76.5583)
    }, action = "append")
}
## https://blog.rstudio.com/2013/06/10/rstudio-cran-mirror/
options(repos = c(CRAN = "https://cloud.r-project.org"))
```

* Es actualizado con bastante frecuencia
* RStudio cheatsheets https://www.rstudio.com/resources/cheatsheets/
  - https://github.com/rstudio/cheatsheets/raw/master/rstudio-ide.pdf
* RStudio projects: usalos para organizar tu código
  - https://github.com/ComunidadBioInfo/cdsb2020/blob/master/presentaciones_flujos-de-trabajo/Trabajando-con-proyectos.pdf



```r
usethis::create_project("~/Desktop/cdsb2021_scRNAseq_notas")
```


```r
## Inicien un archivo para sus notas
usethis::use_r("01-notas.R")
```

O por ejemplo el archivo [01-visualizar-mtcars.R](https://github.com/lcolladotor/rnaseq_2021_notas_en_vivo/blob/master/R/01-visualizar-mtcars.R)


```r
## Creemos el archivo R/01-visualizar-mtcars.R
usethis::use_r("01-visualizar-mtcars.R")
```

con el siguiente contenido:


```r
## Cargar paquetes que usaremos en este código
library("sessioninfo")
library("here")
library("ggplot2")

## Hello world
print("Soy Leo")

## Crear directorio para las figuras
dir.create(here::here("figuras"), showWarnings = FALSE)

## Hacer una imagen de ejemplo
pdf(here::here("figuras", "mtcars_gear_vs_mpg.pdf"),
    useDingbats = FALSE
)
ggplot(mtcars, aes(group = gear, y = mpg)) +
    geom_boxplot()
dev.off()

## Para reproducir mi código
options(width = 120)
sessioninfo::session_info()
```


* Configura `usethis` con GitHub vía https://usethis.r-lib.org/articles/articles/git-credentials.html


```r
## Para poder conectar tu compu con GitHub
usethis::create_github_token() ## Abrirá una página web, escoje un nombre único
## y luego da click en el botón verde al final. Después copia el token
## (son 40 caracteres)

gitcreds::gitcreds_set() ## Ojo, copia el token, no tu password de git!
## Si no, terminaras en la situación descrita en
## https://github.com/r-lib/usethis/issues/1347
```


```r
## Configura tu usuario de GitHub
usethis::edit_git_config()
# [user]
# 	name = Leonardo Collado Torres
# 	email = lcolladotor@gmail.com

## Para inicializar el repositorio de Git
usethis::use_git()

## Para conectar tu repositorio local de Git con los servidores de GitHub
usethis::use_github()
```

Resultado ejemplo: https://github.com/lcolladotor/cdsb2021_scRNAseq_notas. El que hice en vivo está disponible vía https://github.com/lcolladotor/cdsb2021_scRNAseq_notas_en_vivo (o https://github.com/lcolladotor/rnaseq_2021_notas_en_vivo para un ejemplo de febrero 2021).

Una vez que termines, agrega la liga al repositorio con tus notas del curso en el [Google Sheet](https://docs.google.com/spreadsheets/d/13xHCfRb3vATXCFxS1prIA5cYgHNFnzI0GLlcIjtenyw/edit?usp=sharing) del curso. (De ser necesario, pide permisos para editar el archivo.)

## Material del curso

* Pueden descargar la versión estática con `usethis::use_course('ComunidadBioInfo/cdsb2021_scRNAseq')`
* Pueden verlo en línea a través de  [**ComunidadBioInfo.github.io/cdsb2021_scRNAseq**](http://ComunidadBioInfo.github.io/cdsb2021_scRNAseq)
* Pueden **clonarlo** desde GitHub de tal forma que podrán actualizarlo fácilmente usando *git pull*


```bash
git clone https://github.com/ComunidadBioInfo/cdsb2021_scRNAseq.git
## Si tienen su SSH key configurarda pueden usar
## Info sobre SSH keys de GitHub: 
## https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent
git clone git@github.com:ComunidadBioInfo/cdsb2021_scRNAseq.git
```

O desde R con:


```r
## Opción más nueva:
library("gert")
repo <- git_clone(
    "https://github.com/ComunidadBioInfo/cdsb2021_scRNAseq",
    "~/Desktop/cdsb2021_scRNAseq"
)
setwd(repo)

## Otra opción:
git2r::clone(
    "https://github.com/ComunidadBioInfo/cdsb2021_scRNAseq",
    "~/Desktop/cdsb2021_scRNAseq"
)
```

## Detalles de la sesión de R


```r
## Información de la sesión de R
Sys.time()
```

```
## [1] "2021-08-11 04:34:36 UTC"
```

```r
proc.time()
```

```
##    user  system elapsed 
##   0.427   0.144   0.435
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
##  date     2021-08-11                  
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
