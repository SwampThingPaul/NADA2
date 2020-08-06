## Test Environments

* Local: Windows, R version 3.6.1

## R CMD check results (2020-08-04)

### devtools::check() locally:
1 error x | 0 warnings √ | 0 note √

```
> checking examples ... ERROR
  Running examples in 'NADA2-Ex.R' failed
  The error most likely occurred in:
  
  > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
  > ### Name: kenplot
  > ### Title: Plot sensored trend for censored data
  > ### Aliases: kenplot
  > ### Keywords: kendall trend
  > 
  > ### ** Examples
  > 
  > library(NADA) #For example data
  Warning: package 'NADA' was built under R version 3.6.3
  Loading required package: survival
  Warning: package 'survival' was built under R version 3.6.2
  
  Attaching package: 'NADA'
  
  The following object is masked from 'package:stats':
  
      cor
  
  > 
  > # Both y and x are censored
  > data(Golden)
  > with(Golden, kenplot(Blood, BloodCen, Kidney, KidneyCen))
  Error in kenplot(Blood, BloodCen, Kidney, KidneyCen) : 
    could not find function "kenplot"
  Calls: with -> with.default -> eval -> eval
  Execution halted

```

### devtools::check_rhub() locally:
