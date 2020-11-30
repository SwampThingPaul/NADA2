
## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
## Clears Everything...start fresh.
## rm(list=ls(all=T));cat("\014");dev.off()


# Packages used in package ------------------------------------------------
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
library(coin)

# -------------------------------------------------------------------------

## To build the package
library(roxygen2)
library(roxygen2md)
library(pkgdown)

## Process a package
roxygen2::roxygenise()

## Convert from Rd to Markdown
roxygen2md::roxygen2md()


devtools::document()

#devtools::test(); # makes testthat directory and runs tests

devtools::check()

##
devtools::check_rhub()

## Builds package down webpage
pkgdown::build_site()


## Builds package
# devtools::build()



# data(Brumbaugh)
# y.var=Brumbaugh$Hg
# cen.var=Brumbaugh$HgCen
