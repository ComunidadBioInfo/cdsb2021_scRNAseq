## ----install_rsthemes, eval = FALSE-----------------------------------------------------
## remotes::install_github(c(
##     "gadenbuie/rsthemes"
## ))
## remotes::install_cran("suncalc")
## rsthemes::install_rsthemes(include_base16 = TRUE)


## ----r_profile, eval = FALSE------------------------------------------------------------
## usethis::edit_r_profile()
## 
## ## From https://www.garrickadenbuie.com/project/rsthemes/
## if (interactive() && requireNamespace("rsthemes", quietly = TRUE)) {
##     # Set preferred themes if not handled elsewhere..
##     rsthemes::set_theme_light("Solarized Light {rsthemes}") # light theme
##     rsthemes::set_theme_dark("base16 Monokai {rsthemes}") # dark theme
##     rsthemes::set_theme_favorite(c(
##         "Solarized Light {rsthemes}",
##         "base16 Monokai {rsthemes}",
##         "One Dark {rsthemes}"
##     ))
##     # Whenever the R session restarts inside RStudio...
##     setHook("rstudio.sessionInit", function(isNewSession) {
##         # Automatically choose the correct theme based on time of day
##         ## Used rsthemes::geolocate() once
##         rsthemes::use_theme_auto(lat = 39.2891, lon = -76.5583)
##     }, action = "append")
## }
## ## https://blog.rstudio.com/2013/06/10/rstudio-cran-mirror/
## options(repos = c(CRAN = "https://cloud.r-project.org"))


## ----proj, eval = FALSE-----------------------------------------------------------------
## usethis::create_project("~/Desktop/cdsb2021_scRNAseq_notas")


## ----create_setup, eval = FALSE---------------------------------------------------------
## ## Inicien un archivo para sus notas
## usethis::use_r("01-notas.R")


## ----create_01-visualizar-mtcars, eval = FALSE------------------------------------------
## ## Creemos el archivo R/01-visualizar-mtcars.R
## usethis::use_r("01-visualizar-mtcars.R")


## ----vis_mtcars, eval = FALSE-----------------------------------------------------------
## ## Cargar paquetes que usaremos en este código
## library("sessioninfo")
## library("here")
## library("ggplot2")
## 
## ## Hello world
## print("Soy Leo")
## 
## ## Crear directorio para las figuras
## dir.create(here::here("figuras"), showWarnings = FALSE)
## 
## ## Hacer una imagen de ejemplo
## pdf(here::here("figuras", "mtcars_gear_vs_mpg.pdf"),
##     useDingbats = FALSE
## )
## ggplot(mtcars, aes(group = gear, y = mpg)) +
##     geom_boxplot()
## dev.off()
## 
## ## Para reproducir mi código
## options(width = 120)
## sessioninfo::session_info()


## ----use_git_init, eval = FALSE---------------------------------------------------------
## ## Para poder conectar tu compu con GitHub
## usethis::create_github_token() ## Abrirá una página web, escoje un nombre único
## ## y luego da click en el botón verde al final. Después copia el token
## ## (son 40 caracteres)
## 
## gitcreds::gitcreds_set() ## Ojo, copia el token, no tu password de git!
## ## Si no, terminaras en la situación descrita en
## ## https://github.com/r-lib/usethis/issues/1347


## ----use_git_cont, eval = FALSE---------------------------------------------------------
## ## Configura tu usuario de GitHub
## usethis::edit_git_config()
## # [user]
## # 	name = Leonardo Collado Torres
## # 	email = lcolladotor@gmail.com
## 
## ## Para inicializar el repositorio de Git
## usethis::use_git()
## 
## ## Para conectar tu repositorio local de Git con los servidores de GitHub
## usethis::use_github()


## git clone https://github.com/ComunidadBioInfo/cdsb2021_scRNAseq.git

## ## Si tienen su SSH key configurarda pueden usar

## ## Info sobre SSH keys de GitHub:

## ## https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent

## git clone git@github.com:ComunidadBioInfo/cdsb2021_scRNAseq.git


## ----clone_repo, eval = FALSE-----------------------------------------------------------
## ## Opción más nueva:
## library("gert")
## repo <- git_clone(
##     "https://github.com/ComunidadBioInfo/cdsb2021_scRNAseq",
##     "~/Desktop/cdsb2021_scRNAseq"
## )
## setwd(repo)
## 
## ## Otra opción:
## git2r::clone(
##     "https://github.com/ComunidadBioInfo/cdsb2021_scRNAseq",
##     "~/Desktop/cdsb2021_scRNAseq"
## )


## ---------------------------------------------------------------------------------------
## Información de la sesión de R
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

