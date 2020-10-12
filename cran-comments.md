## Test Environments

* Local: Windows, R version 3.6.1

## R CMD check results (2020-08-06)

### devtools::check() locally:
0 errors √ | 0 warnings √ | 0 notes √


### devtools::check_rhub() locally:
0 errors √ | 0 warnings √ | 2 notes x

```
  Non-FOSS package license (file LICENSE)
  
  Possibly mis-spelled words in DESCRIPTION:
    Helsel (13:51)
    Nondetects (14:24)
```


## R CMD check results (2020-10-12)

### devtools::check() locally:
```
> checking R code for possible problems ... NOTE
  cenCompareCdfs: no visible global function definition for 'fitdist'
  cenCompareCdfs: no visible global function definition for 'cdfcomp'
  cenCompareQQ: no visible global function definition for 'distChoose'
  cenQQ: no visible global function definition for 'distChoose'
  Undefined global functions or variables:
    cdfcomp distChoose fitdist
```
0 errors √ | 0 warnings √ | 1 note x
