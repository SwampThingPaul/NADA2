## Internal Data
## From R-package NADA Atra dataset
## Atrazine concentrations

atrazine=read.csv("./data-raw/NADA_atrazine.csv")
usethis::use_data(atrazine,internal=F,overwrite=T)
