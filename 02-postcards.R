## ---- eval=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------
## # se descargó previamente, así que solo se carga
## library("here") # busca la raiz del proyecto en el que se encuentre


## ---- eval=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------
## here::here()


## ---- eval=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------
## getwd() # regresa la path en donde nos encontramos
## setwd("direccion/deseada") # nos lleva a la path indicada


## ---- eval=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------
## getwd() # para checar en donde nos encontramos
## here::here() # para checar dónde te encuentras
## 
## # nos movemos al subdirectorio R
## setwd(here::here("R")) # podemos cambiar de directorio, aun así `here está en la raíz


## ---- eval=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------
## # como ejemplo: vamos a guardar datos en archivo y cargarlos
## a <- 1
## c <- 23
## save(a, c, file = here::here("datos-prueba.RData"))
## # rm(a,c)
## load(here::here("datos-prueba.RData"))
## 
## # creamos un directorio
## dir.create(here::here("subdirectorio"), showWarnings = FALSE)
## # podemos crear un archivo, indicando el subdirectorio, (en este caso el primer argumento)
## file.create(here::here("subdirectorio", "nombrearchivo"))
## # abrimos el nuevo archivo creado
## file.show(here::here("subdirectorio", "nombrearchivo")) # podemos editarlo!!
## 
## # por ejemplo si quisieramos ver nuestros archivos del directorio
## list.files(here::here(), recursive = TRUE)


## ---- eval=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------
## usethis::use_r("notas-prueba.R") # no importando en qué path estemos


## ---- eval=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------
## # paquetes que vamos a requerir
## install.packages(c("gitcreds", "gert", "gh"))
## # cargarlos de manera separada
## library("gitcreds")
## library("gert")
## library("gh")


## ---- eval=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------
## # Para iniciar conexión con GitHub
## usethis::create_github_token() # redirige a github donde eligiras nombre especifico del token
## # copia el token para después ingresarlo con gitcreds_set()
## gitcreds::gitcreds_set() # aquí colocas el token (NO tu contraseña de github!!!)


## ---- eval=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------
## # Configurar usuario de gitHub
## usethis::edit_git_config() # que abre el archivo .gitconfig
## # colocaremos nombre y correo de cuenta de github. SOLO borrar los # y respetar los demas espacios
## # [user]
## #   name = N O M B R E
## #   email = correodeGithub


## ---- eval=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------
## # inicializar el repositorio de Git
## usethis::use_git() #
## 
## # conectar tu repositorio local de Git con los servidores de GitHub
## usethis::use_github()


## ---- eval=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------
## gh::gh_whoami() # para checar cómo quedó la configuración


## ---- eval=FALSE, warning=FALSE, message=FALSE-------------------------------------------------------
## # escribimos un nuevo archivo, volvemos a usar here::here para especificar path
## writeLines("hola", here::here("R", "prueba-here.R"))
## # otra manera es usar use_r
## usethis::use_r("archivo-prueba-github.R") # añade archivo al directorio R del proyecto actual
## 
## # Por ejemplo podríamos probar añadir algo nuevo
## gert::git_add("R/archivo-prueba-github.R")
## 
## # añadimos commit de lo que se hizo
## gert::git_commit("se subio archivo prueba")
## 
## # nos da info de los commits
## gert::git_log()
## 
## # sube tus cambios del repo local a los de github
## gert::git_push() # COMANDO IMPORTANTE


## ----postcards_proj, eval = FALSE--------------------------------------------------------------------
## ## Creen el RStudio project. Es MUY importante que el usuario debe sea igual que en github
## usethis::create_project("Su_Usuario.github.io")


## ----echo=FALSE, fig.cap="git user", out.width = "40%"-----------------------------------------------
knitr::include_graphics("img_postcard/proyecto_github_user.png")


## ----postcards_git, eval = FALSE---------------------------------------------------------------------
## ## Configura Git y GitHub
## 
## # Con use_git() preguntará si desean hacer un commit, y después pedirá reiniciar Rstudio para que obtengan un nuevo botón llamado "git()"
## usethis::use_git()


## ----echo=FALSE, fig.cap="button_git", out.width = "40%"---------------------------------------------
knitr::include_graphics("img_postcard/activa_boton_git.png")


## ----postcards_github, eval = FALSE------------------------------------------------------------------
## usethis::use_github()


## ----postcards_create, eval = FALSE------------------------------------------------------------------
## ## Solo uno de estos, de acuerdo al templado que más les gustó
## postcards::create_postcard(template = "jolla")
## postcards::create_postcard(template = "jolla-blue")
## postcards::create_postcard(template = "trestles")
## postcards::create_postcard(template = "onofre")


## ----postcards_git2, eval = FALSE--------------------------------------------------------------------
## ## Configura Git y GitHub
## 
## # Con use_git() preguntará si desean hacer un commit, y después pedirá reiniciar Rstudio para que obtengan un nuevo botón llamado "git()"
## usethis::use_git()


## ----echo=FALSE, fig.cap="button_git", out.width = "40%"---------------------------------------------
knitr::include_graphics("img_postcard/activa_boton_git.png")


## ----echo=FALSE, fig.cap="button_branch", out.width = "40%"------------------------------------------
knitr::include_graphics("img_postcard/branch.png")


## ----echo=FALSE, fig.cap="add_file", out.width = "40%"-----------------------------------------------
knitr::include_graphics("img_postcard/add_stag.png")


## ----echo=FALSE, fig.cap="button_git", out.width = "80%"---------------------------------------------
knitr::include_graphics("img_postcard/commit_message.png")


## ----postcards_pull, eval = FALSE--------------------------------------------------------------------
## ## Guardamos los archivos nuevos con el commit
## gert::git_commit_all("mensaje sobre el commit")
## 
## ## Subimos los archivos a github
## gert::git_push


## ----fawesome_yaml, eval = FALSE---------------------------------------------------------------------
## 
## 
## # Utlilizando `r fontawesome::fa("font-awesome-logo-full", fill = "forestgreen")` en código YAML
## 


## ----echo=FALSE, fig.cap="y_fawesome", out.width = "40%"---------------------------------------------
knitr::include_graphics("img_postcard/yaml_awesome.png")


## ----echo=FALSE, fig.cap="fawesome", out.width = "40%"-----------------------------------------------
knitr::include_graphics("img_postcard/awesome.png")


## ----fawesome_text, eval = FALSE---------------------------------------------------------------------
## 
## 
## # Utlilizando `r fontawesome::fa("font-awesome-logo-full", fill = "forestgreen")` en el texto
## 


## ----echo=FALSE, fig.cap="tfawesome", out.width = "40%"----------------------------------------------
knitr::include_graphics("img_postcard/text_awesome.png")


## ----------------------------------------------------------------------------------------------------
## Información de la sesión de R
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

