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
  ROSci: no visible global function definition for 'qt'
  ROSci: no visible global function definition for 'sd'
  anosimPlot: no visible global function definition for 'hist'
  binaryClust: no visible global function definition for 'hclust'
  binaryClust: no visible global function definition for 'cutree'
  binaryClust: no visible global function definition for 'rect.hclust'
  binaryMDS: no visible global function definition for 'ordiplot'
  binaryMDS: no visible global function definition for 'points'
  cboxplot: no visible global function definition for 'adjustcolor'
  cboxplot: no visible global function definition for 'ros'
  cboxplot: no visible global function definition for 'polygon'
  cboxplot: no visible global function definition for 'sd'
  cen1way: no visible global function definition for 'pchisq'
  cen2means: no visible global function definition for 'pchisq'
  cen2means: no visible global function definition for 'predict'
  cenTolInt: no visible global function definition for 'eqlnormCensored'
  cenTolInt: no visible global function definition for 'eqnormCensored'
  cen_paired: no visible global function definition for 'plotdistcens'
  cen_signedranktest: no visible global function definition for
    'wilcoxsign_test'
  cen_signedranktest: no visible global function definition for
    'statistic'
  cen_signedranktest: no visible global function definition for 'pvalue'
  cen_signtest: no visible global function definition for 'binom.test'
  cen_signtest: no visible global function definition for 'pbinom'
  cen_signtest: no visible global function definition for 'dbinom'
  cenanova: no visible global function definition for 'pchisq'
  cenanova: no visible global function definition for 'mcp'
  cenanova: no visible global function definition for 'predict'
  cenanova: no visible global function definition for 'residuals'
  cenpermanova: no visible global function definition for 'cenros'
  censeaken: no visible global function definition for 'hist'
  cfit: no visible global function definition for 'qt'
  cfit: no visible global function definition for 'quantile'
  plotcdf: no visible global function definition for 'box'
  plotcdf: no visible global function definition for 'axis'
  ppw.test : PPW.test: no visible binding for global variable
    'na.exclude'
  ppw.test: no visible global function definition for 'pnorm'
  uMDS: no visible global function definition for 'dist'
  uMDS: no visible global function definition for 'ordiplot'
  uMDS: no visible global function definition for 'points'
  Undefined global functions or variables:
    adjustcolor axis binom.test box cenros cutree dbinom dist
    eqlnormCensored eqnormCensored hclust hist mcp na.exclude ordiplot
    pbinom pchisq plotdistcens pnorm points polygon predict pvalue qt
    quantile rect.hclust residuals ros sd statistic wilcoxsign_test
  Consider adding
    importFrom("grDevices", "adjustcolor")
    importFrom("graphics", "axis", "box", "hist", "points", "polygon")
    importFrom("stats", "binom.test", "cutree", "dbinom", "dist", "hclust",
               "na.exclude", "pbinom", "pchisq", "pnorm", "predict", "qt",
               "quantile", "rect.hclust", "residuals", "sd")
  to your NAMESPACE file.

```
