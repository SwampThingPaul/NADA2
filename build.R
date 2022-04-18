
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
# https://johnmuschelli.com/smi_2019/index.html#78
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

# NEWs file
# usethis::use_news_md()

## Builds package down webpage
pkgdown::build_site()

## Builds package
# devtools::build()


# -------------------------------------------------------------------------

library(NADA2)


# centrend

y.var=Brumbaugh$Hg
y.cens=Brumbaugh$HgCen
x.var=Brumbaugh$SedTotHg
time.var=1:nrow(Brumbaugh)
Smooth = "cs"
link = "identity"



####
test <- Brumbaugh

# Estimate the mean with ROS
Erg <- NADA::ros(test$Hg, test$HgCen)
summary(Erg)

mean(Erg); # produces NAs
mean(Erg$modeled); # value of 0.3555983

Erg2 <- NADA::cenros(test$Hg, test$HgCen)
summary(Erg2)

mean(Erg); # produces NAs
mean(Erg2$modeled); # value of 0.3555983

# Estimate confidence interval of mean
ROSci(Erg)
