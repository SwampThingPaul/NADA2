## Test Environments

* Local: Windows, R version 3.6.1

## R CMD check results (2020-08-06)

### devtools::check() locally:
1 error x | 0 warnings √ | 0 note √

```
> checking examples ... ERROR
  Running examples in 'NADA2-Ex.R' failed
  The error most likely occurred in:
  
  > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
  > ### Name: plotcdf
  > ### Title: Censored Empirical Cumulative Distribution Function
  > ### Aliases: plotcdf
  > ### Keywords: CDF
  > 
  > ### ** Examples
  > 
  > data(PbHeron)
  > 
  > # with groups
  > with(PbHeron,plotcdf(Liver,LiverCen,DosageGroup))
  Error in as.double(y) : 
    cannot coerce type 'S4' to vector of type 'double'
  Calls: with ... plotcdf -> stripaxes -> plot -> plot.default -> xy.coords
  Execution halted


```

### devtools::check_rhub() locally:
