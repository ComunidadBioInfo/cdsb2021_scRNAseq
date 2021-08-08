# Ejercicio usando `usethis`, `here` y `postcards`

## Ejercicio postcards
* Similar a https://pages.github.com/
* `postcards` tiene 4 templados de páginas web https://github.com/seankross/postcards
* Tu página web debe describir decir algo sobre ti, tus intereses, y tus proyectos además de cómo contactarte
  - Ejemplo https://amy-peterson.github.io/ vía https://github.com/amy-peterson/amy-peterson.github.com
  - http://jtleek.com/ vía https://github.com/jtleek/jtleek.github.io
  - http://aejaffe.com/ vía https://github.com/andrewejaffe/andrewejaffe.github.io


### Crear el repositorio
* Podemos crearlo desde RStudio o desde github.com

**(opción 1) Desde RStudio**

<div class="alert alert-block alert-danger">
<b>¡Cuidado!:</b> Antes de crear un proyecto, revisen dónde están parados (`getwd()`) en su directorio y dónde quieren que se cree
</div>


```r
## Creen el RStudio project. Es MUY importante que el usuario debe sea igual que en github
usethis::create_project("Su_Usuario.github.io")
```

Nuevo proyecto :
<img src="img_postcard/proyecto_github_user.png" alt="git_user" width="300"/>

```r
## Configura Git y GitHub

# Con use_git() preguntará si desean hacer un commit, y después pedirá reiniciar Rstudio para que obtengan un nuevo botón llamado "git()"
usethis::use_git()
```

Nuevo botón: <img src="img_postcard/activa_boton_git.png" alt="button_git" width="300"/>



```r
usethis::use_github()
```


Creen su templado usando `postcards`. Va a crear un archivo `index.Rmd`


```r
## Solo uno de estos, de acuerdo al templado que más les gustó
postcards::create_postcard(template = "jolla")
postcards::create_postcard(template = "jolla-blue")
postcards::create_postcard(template = "trestles")
postcards::create_postcard(template = "onofre")
```

**(opción 2) Desde github**

* Creen un nuevo repositorio, público y **sin** archivo README en https://github.com/new llamado "usuario.github.io" con su **nombre exacto** en github

<div class="alert alert-warning">
  <strong>¡Cuidado!</strong> El repositorio debe ser **público** y sin README
</div>

* Creen un nuevo proyecto en RStudio: 
  File > New_project > New directory > Postcards Website 
  
* Elijan el templado que más les gustó

* Ya con el proyecto creado, hay que configurar git y github


```r
## Configura Git y GitHub

# Con use_git() preguntará si desean hacer un commit, y después pedirá reiniciar Rstudio para que obtengan un nuevo botón llamado "git()"
usethis::use_git()
```

  
  
Nuevo botón: <img src="img_postcard/activa_boton_git.png" alt="button_git" width="300"/>

* Ahora que tienen el botón Git, hagan click y en la esquina derecha habrá un símbolo con dos  rectángulos morados y un rombo blanco, denle click.


 <img src="img_postcard/branch.png" alt="button_brach" width="300"/>

* Ahora el botón *Add Remote* y ahí podrán nombrar este acceso remoto como gusten, y agregar la URL de su repositorio en github.

* Da click en *Add* y después asignen el nombre de rama **master**

* Ahora pueden crear la rama, y sobreescribir el acceso cuando se los pregunte.

-----------------------------


-----------------------------

### Modificar y subir a github nuestro postcard

Ya que hayan creado con cualquiera de las 2 opciones anteriores pueden continuar:

* Llenen su información usando el formato `Markdown`. Por ejemplo https://github.com/andrewejaffe/andrewejaffe.github.io/blob/master/index.Rmd#L17-L31.
* Agreguen sus perfiles estilo https://github.com/andrewejaffe/andrewejaffe.github.io/blob/master/index.Rmd#L7-L12
* Den click en el botón azul de `knit` en RStudio. Es equivalente a `rmarkdown::render("index.Rmd")`. Esto creará el archivo `index.html`.

* Hagan un `commit` para guardar los archivos nuevos incluyendo `index.html` y luego un `push` para subir los archivos a GitHub con alguna de las siguientes dos maneras:

**(opción 1) Botón de Git**
* Para guardar los archivos nuevos, incluyendo `index.html`, debemos hacer un `commit`. Podemos hacerlo con el nuevo botón de git, primero seleccionando los archivos: 
<img src="img_postcard/add_stag.png" alt="add_file" width="300"/>

Cuando hayamos seleccionado **todos** los archivos, veremos que la columna Status cambia a una "A" de agregado o added y podemos darle al botón `Commit` justo arriba de Status. Esto abrirá una nueva pestaña donde podremos poner un mensaje sobre nuestro `commit` y después darle al botón `Commit`.


<img src="img_postcard/commit_message.png" alt="add_file" width="800"  class="center"/>

Una vez terminado, en esa misma pantalla podemos darle un `push` para subir los archivos a GitHub con el botón de  `Push` con una flecha verde arriba de *Commit message*.


**(opcion 2) Línea de comandos**
Otra manera de hacer es vía línea de comandos, primero pueden agregar los archivos con `gert::git_add()`  o hacer directamente un commit de todos los archivos y luego un pull:


```r
## Guardamos los archivos nuevos con el commit
gert::git_commit_all('mensaje sobre el commit')

## Subimos los archivos a github
gert::git_push
```

* <span style="color:DodgerBlue">**(extra)** </span>. Pueden copiar y pegar emojis en sus páginas o utilizar [fontawesome](https://github.com/rstudio/fontawesome) para agregar diferentes símbolos (como github o twitter):

  - En código YAML

```r
# Utlilizando `r fontawesome::fa("font-awesome-logo-full", fill = "forestgreen")` en código YAML
```


  
  <img src="img_postcard/yaml_awesome.png" alt="y_fawesome" width="800"/>
  
  
  
  Se ve así:
  
  
  <img src="img_postcard/awesome.png" alt="fawesome" width="500"/>

* 
   - En el texto



```r
# Utlilizando `r fontawesome::fa("font-awesome-logo-full", fill = "forestgreen")` en el texto
```

<img src="img_postcard/text_awesome.png" alt="tfawesome" width="300"/>



* <span style="color:DodgerBlue">**(opcional)**</span>. Anuncien su nueva página web en Twitter usando el hashtag `#rstats` y/o etiquen al autor de `postcards` https://twitter.com/seankross. Pueden después incluir su página web en su introducción en el canal `#bienvenida` del Slack de la CDSB ^^.

## Detalles de la sesión de R


```r
## Información de la sesión de R
Sys.time()
```

```
## [1] "2021-08-08 21:00:25 UTC"
```

```r
proc.time()
```

```
##    user  system elapsed 
##   0.438   0.139   0.449
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
##  date     2021-08-08                  
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
##  rmarkdown     2.9     2021-06-15 [1] RSPM (R 4.1.0)
##  sass          0.4.0   2021-05-12 [1] RSPM (R 4.1.0)
##  sessioninfo   1.1.1   2018-11-05 [2] RSPM (R 4.1.0)
##  stringi       1.7.3   2021-07-16 [2] RSPM (R 4.1.0)
##  stringr       1.4.0   2019-02-10 [2] RSPM (R 4.1.0)
##  withr         2.4.2   2021-04-18 [2] RSPM (R 4.1.0)
##  xfun          0.24    2021-06-15 [2] RSPM (R 4.1.0)
## 
## [1] /__w/_temp/Library
## [2] /usr/local/lib/R/site-library
## [3] /usr/local/lib/R/library
```

## Patrocinadores {-}

Agradecemos a nuestros patrocinadores:

<a href="https://comunidadbioinfo.github.io/es/post/cs_and_s_event_fund_award/#.YJH-wbVKj8A"><img src="https://comunidadbioinfo.github.io/post/2021-01-27-cs_and_s_event_fund_award/spanish_cs_and_s_award.jpeg" width="400px" align="center"/></a>

<a href="https://www.r-consortium.org/"><img src="https://www.r-consortium.org/wp-content/uploads/sites/13/2016/09/RConsortium_Horizontal_Pantone.png" width="400px" align="center"/></a>
