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
  > ### Name: cenCompareQQ
  > ### Title: Censored Q-Q Plot comparison
  > ### Aliases: cenCompareQQ
  > 
  > ### ** Examples
  > 
  > data(Brumbaugh)
  > 
  > cenCompareQQ(Brumbaugh$Hg,Brumbaugh$HgCen)
  Error in as.environment(where) : 
    no item called "package:EnvStats" on the search list
  Calls: cenCompareQQ ... envstatsDistChooseCensored -> gofTestCensored -> exists
  Execution halted

```
