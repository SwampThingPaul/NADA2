## Test Environments

* Local: Windows, R version 3.6.1

## R CMD check results (2020-07-23)

### devtools::check() locally:
1 error x | 0 warnings √ | 1 note x

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
  > cenCompareQQ(Brumbaugh$Hg,Brumbaugh$HgCen)
  Error in as.environment(where) : 
    no item called "package:EnvStats" on the search list
  Calls: cenCompareQQ ... envstatsDistChooseCensored -> gofTestCensored -> exists
  Execution halted

> checking R code for possible problems ... NOTE

  ppw.test : PPW.test: no visible binding for global variable
    'na.exclude'
  ppw.test: no visible global function definition for 'pnorm'
  uMDS: no visible global function definition for 'dist'
  Undefined global functions or variables:
    axis binom.test box cenros dbinom dist mcp na.exclude pbinom
    plotdistcens pnorm pvalue quantile residuals statistic
    wilcoxsign_test
  Consider adding
    importFrom("graphics", "axis", "box")
    importFrom("stats", "binom.test", "dbinom", "dist", "na.exclude",
               "pbinom", "pnorm", "quantile", "residuals")
  to your NAMESPACE file.


```
