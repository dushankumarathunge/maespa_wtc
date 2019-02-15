# do this once
# devtools::install_bitbucket("remkoduursma/plotby")
library(plotBy)

r <- require(HIEv)
if(!r)stop("Install the HIEv R package from bitbucket.org/remkoduursma/hiev")





#---------------------------------------------------------------------
#- function to load a package, and install it if necessary
Library <- function(pkg, ...){
  
  PACK <- .packages(all.available=TRUE)
  pkgc <- deparse(substitute(pkg))
  
  if(pkgc %in% PACK){
    library(pkgc, character.only=TRUE)
  } else {
    install.packages(pkgc, ...)
    library(pkgc, character.only=TRUE)
  }
  
}
#---------------------------------------------------------------------





#---------------------------------------------------------------------
#- load all the libraries (and install them if needed)
Library(doBy)
Library(magicaxis)
Library(RColorBrewer)
Library(propagate)
Library(gplots)
Library(readxl)
Library(maps)
Library(mapdata)
Library(rgeos)
Library(sp)
Library(raster)
Library(nlstools)
Library(rgl)
Library(mgcv)
Library(scales)
Library(data.table)
Library(dplyr)
Library(dismo)
Library(plantecophys)
Library(plotrix)
Library(kableExtra)
Library(Maeswrap)
Library(HIEv)
Library(gplots)
Library(zoo)
Library(doBy)
# Library(plover)
#---------------------------------------------------------------------






#---------------------------------------------------------------------
#- load the custom functions that do most of the heavy lifting

# source("R/generic_functions.R")
# source("R/GREAT_functions.R")
#source("R/gamplotfunctions.R")

#---------------------------------------------------------------------
