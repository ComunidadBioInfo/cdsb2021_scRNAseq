## ----install, eval = FALSE--------------------------------------------------------------------------
## ## Para instalar paquetes
## if (!requireNamespace("remotes", quietly = TRUE)) {
##     install.packages("remotes")
## }
## 
## ## Para instalar paquetes de Bioconductor
## remotes::install_cran("BiocManager")
## BiocManager::version()
## # El anterior comando debe mostrar que estás usando la versión 3.13
## 
## ## Instala los paquetes de R que necesitamos


## ----auto_deps, echo = FALSE------------------------------------------------------------------------
cat(paste0(
    'BiocManager::install(c("',
    paste(
        sort(remotes::local_package_deps()),
        collapse = '", "'
    ),
    '"))'
))


## ----echo=FALSE, fig.cap="Al crear un nuevo proyecto, seleccionen la opción de _Version Control_ (la tercera).", out.width="60%"----
knitr::include_graphics("img/clone_version_control.png")


## ----echo=FALSE, fig.cap="Selecciona la opción de `Git` (la primera).", out.width="60%"-------------
knitr::include_graphics("img/clone_choose_git.png")


## ----echo=FALSE, fig.cap="Especifica que el _Repository URL_ es `https://github.com/ComunidadBioInfo/cdsb2021_scRNAseq.git`.", out.width="60%"----
knitr::include_graphics("img/clone_add_info.png")


## ----session_info-----------------------------------------------------------------------------------
options(width = 120)
pkgs <- installed.packages()[, "Package"]
sessioninfo::session_info(pkgs, include_base = TRUE)

