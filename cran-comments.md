## Test Environments

* Local: Windows, R version 3.6.1

## R CMD check results (2020-07-23)

### devtools::check() locally:
1 errors | 3 warnings | 1 notes

```

> checking whether package 'NADA2' can be installed ... WARNING
  See below...

> checking dependencies in R code ... WARNING
  Package in Depends field not imported from: 'knitr'
    These packages need to be imported from (in the NAMESPACE file)
    for when this namespace is loaded but not attached.
  Unexported object imported by a ':::' call: 'survival:::survmean'
    See the note in ?`:::` about the use of this operator.
    Including base/recommended package(s):
    'survival'

> checking Rd \usage sections ... WARNING
  Undocumented arguments in documentation object 'equivalent_n'
    'y.cen'
  Documented arguments not in \usage in documentation object 'equivalent_n':
    'cen.var'
  
  Functions with \usage entries need to have the appropriate \alias
  entries, and all their arguments documented.
  The \usage entries must correspond to syntactically valid R code.
  See chapter 'Writing R documentation files' in the 'Writing R
  Extensions' manual.

> checking R code for possible problems ... NOTE
  ATS: no visible global function definition for 'na.omit'
  ROSci: no visible global function definition for 'qt'
  ROSci: no visible global function definition for 'sd'
  Usc: no visible global function definition for 'na.omit'
  anosimPlot: no visible global function definition for 'hist'
  binaryClust: no visible global function definition for 'hclust'
  binaryClust: no visible global function definition for 'cutree'
  binaryClust: no visible global function definition for 'rect.hclust'
  binaryMDS: no visible global function definition for 'ordiplot'
  binaryMDS: no visible global function definition for 'points'
  cboxplot: no visible global function definition for 'adjustcolor'
  cboxplot: no visible global function definition for 'na.omit'
  cboxplot: no visible global function definition for 'ros'
  cboxplot: no visible global function definition for 'polygon'
  cboxplot: no visible global function definition for 'sd'
  cen1way: no visible global function definition for 'pchisq'
  cen1way: no visible global function definition for 'pairwise_survdiff'
  cen2means: no visible global function definition for 'pchisq'
  cen2means: no visible global function definition for 'predict'
  cenTolInt: no visible global function definition for 'eqlnormCensored'
  cenTolInt: no visible global function definition for 'eqnormCensored'
  cen_paired: no visible global function definition for 'na.omit'
  cen_paired: no visible global function definition for 'plotdistcens'
  cen_signedranktest: no visible global function definition for 'na.omit'
  cen_signedranktest: no visible global function definition for
    'wilcoxsign_test'
  cen_signedranktest: no visible global function definition for
    'statistic'
  cen_signedranktest: no visible global function definition for 'pvalue'
  cen_signtest: no visible global function definition for 'na.omit'
  cen_signtest: no visible global function definition for 'binom.test'
  cen_signtest: no visible global function definition for 'pbinom'
  cen_signtest: no visible global function definition for 'dbinom'
  cenanova: no visible global function definition for 'pchisq'
  cenanova: no visible global function definition for 'mcp'
  cenanova: no visible global function definition for 'predict'
  cenanova: no visible global function definition for 'residuals'
  cencorreg: no visible global function definition for 'na.omit'
  cenperm2: no visible global function definition for 'na.omit'
  cenpermanova: no visible global function definition for 'na.omit'
  cenpermanova: no visible global function definition for 'cenros'
  censeaken: no visible global function definition for 'na.omit'
  censeaken: no visible global function definition for 'median'
  censeaken: no visible global function definition for 'hist'
  centrend: no visible global function definition for 'na.omit'
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
  uscore: no visible global function definition for 'na.omit'
  uscores: no visible global function definition for 'na.omit'
  uscoresi: no visible global function definition for 'na.omit'
  Undefined global functions or variables:
    adjustcolor axis binom.test box cenros cutree dbinom dist
    eqlnormCensored eqnormCensored hclust hist mcp median na.exclude
    na.omit ordiplot pairwise_survdiff pbinom pchisq plotdistcens pnorm
    points polygon predict pvalue qt quantile rect.hclust residuals ros
    sd statistic wilcoxsign_test
  Consider adding
    importFrom("grDevices", "adjustcolor")
    importFrom("graphics", "axis", "box", "hist", "points", "polygon")
    importFrom("stats", "binom.test", "cutree", "dbinom", "dist", "hclust",
               "median", "na.exclude", "na.omit", "pbinom", "pchisq",
               "pnorm", "predict", "qt", "quantile", "rect.hclust",
               "residuals", "sd")
  to your NAMESPACE file.
```
