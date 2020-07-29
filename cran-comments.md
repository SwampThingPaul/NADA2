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
  > ### Name: cenPredInt
  > ### Title: Prediction interval for censored data
  > ### Aliases: cenPredInt
  > ### Keywords: interval prediction
  > 
  > ### ** Examples
  > 
  > 
  > library(NADA) #for example data
  Warning: package 'NADA' was built under R version 3.6.3
  Loading required package: survival
  Warning: package 'survival' was built under R version 3.6.2
  
  Attaching package: 'NADA'
  
  The following object is masked from 'package:stats':
  
      cor
  
  > data(Golden)
  > 
  > # Default
  > cenPredInt(Golden$Liver,Golden$LiverCen)
  95% Prediction Limits
    Distribution       95% LPL  95% UPL
  1    Lognormal   0.001464239 28.13403
  2        Gamma   0.000000000 24.93266
  3       Normal -24.332274410 30.61011
  > 
  > # User defined confidence coefficient
  > cenPredInt(Golden$Liver,Golden$LiverCen, conf=0.5)
  50% Prediction Limits
    Distribution       50% LPL   50% UPL
  1    Lognormal  0.0393247897  1.047557
  2        Gamma  0.0001337408  3.271342
  3       Normal -6.0029976395 12.280833
  > 
  > # User defined confidence coefficient outside of acceptable range
  > cenPredInt(Golden$Liver,Golden$LiverCen, conf=1.1)
  Error in cenPredInt(Golden$Liver, Golden$LiverCen, conf = 1.1) : 
    Please select a confidence coefficient between 0 and 1.
  Execution halted

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
    eqlnormCensored eqnormCensored hclust hist mcp na.exclude na.omit
    ordiplot pbinom pchisq plotdistcens pnorm points polygon predict
    pvalue qt quantile rect.hclust residuals ros sd statistic
    wilcoxsign_test
  Consider adding
    importFrom("grDevices", "adjustcolor")
    importFrom("graphics", "axis", "box", "hist", "points", "polygon")
    importFrom("stats", "binom.test", "cutree", "dbinom", "dist", "hclust",
               "na.exclude", "na.omit", "pbinom", "pchisq", "pnorm",
               "predict", "qt", "quantile", "rect.hclust", "residuals",
               "sd")
  to your NAMESPACE file.


```
