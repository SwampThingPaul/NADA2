## Internal Data
## Thiamethoxam concentrations in pollen from the Ontario Pollen Monitoring
## Network (https://data.ontario.ca/en/dataset/pollen-monitoring-network-study).
## Downloaded in 2019.

PollenThia=load("./data-raw/PollenThia.RData")
usethis::use_data(PollenThia,internal=F,overwrite=T)
