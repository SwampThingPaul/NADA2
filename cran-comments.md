## Test Environments

* Local: Windows, R version 3.6.1

## R CMD check results (2020-07-23)

### devtools::check() locally:
1 errors | 0 warnings | 1 notes

```

> checking examples ... ERROR
  Running examples in 'NADA2-Ex.R' failed
  The error most likely occurred in:
  
  > base::assign(".ptime", proc.time(), pos = "CheckExEnv")
  > ### Name: binaryMDS
  > ### Title: Plot Nonmetric Multidimensional Scaling of censored data
  > ### Aliases: binaryMDS
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
  
  > data(Golden)
  > 
  > # without group specified
  > binaryMDS(Golden[,4:15])
  Run 0 stress 9.969016e-05 
  Run 1 stress 7.055086e-05 
  ... New best solution
  ... Procrustes: rmse 0.1002823  max resid 0.240009 
  Run 2 stress 8.26136e-05 
  ... Procrustes: rmse 0.06052965  max resid 0.1754196 
  Run 3 stress 0 
  ... New best solution
  ... Procrustes: rmse 0.04647656  max resid 0.1273565 
  Run 4 stress 7.145941e-05 
  ... Procrustes: rmse 0.1191691  max resid 0.3438194 
  Run 5 stress 0.0006572667 
  Run 6 stress 0 
  ... Procrustes: rmse 0.1319002  max resid 0.3643102 
  Run 7 stress 4.878085e-06 
  ... Procrustes: rmse 0.1223374  max resid 0.3510204 
  Run 8 stress 9.784276e-05 
  ... Procrustes: rmse 0.1162553  max resid 0.3477434 
  Run 9 stress 9.913202e-05 
  ... Procrustes: rmse 0.1128221  max resid 0.1837219 
  Run 10 stress 0.0579766 
  Run 11 stress 9.074527e-05 
  ... Procrustes: rmse 0.1359424  max resid 0.3639362 
  Run 12 stress 0.0002766434 
  ... Procrustes: rmse 0.07981133  max resid 0.1574113 
  Run 13 stress 9.626594e-05 
  ... Procrustes: rmse 0.1168768  max resid 0.3486435 
  Run 14 stress 0.002717952 
  Run 15 stress 5.69538e-05 
  ... Procrustes: rmse 0.1340435  max resid 0.3538703 
  Run 16 stress 0.05797344 
  Run 17 stress 7.44569e-05 
  ... Procrustes: rmse 0.131565  max resid 0.3692676 
  Run 18 stress 0.001999098 
  Run 19 stress 5.414727e-05 
  ... Procrustes: rmse 0.1152273  max resid 0.3485584 
  Run 20 stress 6.77439e-05 
  ... Procrustes: rmse 0.05314647  max resid 0.1178893 

> checking R code for possible problems ... NOTE
  ATS: no visible global function definition for 'na.omit'
  NADA2.survmean: no visible global function definition for 'median'
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
