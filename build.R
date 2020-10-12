
library(EnvStats)
library(fitdistrplus)
library(Kendall)
library(multcomp)
library(NADA)
library(perm)
library(survminer)
library(mgcv)
library(cenGAM)
library(vegan)


## Double check setMethod and class
# http://r-pkgs.had.co.nz/man.html#man-classes


## To build the package
library(roxygen2)
library(roxygen2md)
library(pkgdown)

## Process a package
roxygen2::roxygenise()

## Convert from Rd to Markdown
roxygen2md::roxygen2md()

devtools::document()

# devtools::test(); # makes testthat directory and runs tests

devtools::check()

##
devtools::check_rhub()

## Builds package down webpage
pkgdown::build_site()


## Builds package
# devtools::build()



data(Brumbaugh)
y.var=Brumbaugh$Hg
cen.var=Brumbaugh$HgCen
