## ----install, eval = FALSE----------------------------------------------------------------------------
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


## ----auto_deps, echo = FALSE--------------------------------------------------------------------------
cat(paste0(
    'BiocManager::install(c("',
    paste(
        sort(remotes::local_package_deps()),
        collapse = '", "'
    ),
    '"))'
))


## ----session_info-------------------------------------------------------------------------------------
options(width = 120)
pkgs <- installed.packages()[, "Package"]
sessioninfo::session_info(pkgs, include_base = TRUE)

