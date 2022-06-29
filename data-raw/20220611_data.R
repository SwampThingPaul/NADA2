## New Internal Data - Internal Data

## Moved to Data folder
# load("./data-raw/Example1.RData")
# usethis::use_data(Example1,internal=F,overwrite=T)
#
# load("./data-raw/Example2.RData")
# usethis::use_data(Example2,internal=F,overwrite=T)
#
# load("./data-raw/Example3.RData")
# usethis::use_data(Example3,internal=F,overwrite=T)

load("./data-raw/Markers.RData")
usethis::use_data(Markers,internal=F,overwrite=T)

load("./data-raw/ReconLogistic.RData")
usethis::use_data(ReconLogistic,internal=F,overwrite=T)

load("./data-raw/TCE2.RData")
usethis::use_data(TCE2,internal=F,overwrite=T)

load("./data-raw/Gales_Creek.RData")
colnames(Gales_Creek)=c("Date","Yr","Month","Season","Descript","Units","Flag","TCr","discharge","CrND","dectime")
usethis::use_data(Gales_Creek,internal=F,overwrite=T)
