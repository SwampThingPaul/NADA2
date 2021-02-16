
# NADA2 <img src="man/figures/logo.png" align="right" alt="" width="120" />

## Table of Contents

-   [Introduction](#intro)
-   [Citing package](#cite)
-   [Installation](#install)

## Introduction <a name="intro"></a>

For many environmental professionals, the way to deal with “nondetects”
is less than obvious. Values below detection or reporting limits result
from measuring trace amounts of a variety of organic and inorganic
chemicals. This package includes methods for plotting, computing summary
statistics, hypothesis tests, and regression models for data with one or
more detection limits. Most are adaptations of or direct use of existing
methods found in survival analysis (such as in the survival package) and
in cancer research (such as in the interval package).

This package is a series of new procedures or improvements on current
ones in the `NADA` package developed through a series of webinars,
online courses, and updates to Helsel (2004) and Helsel (2011).

-   Helsel, D.R., 2004. Nondetects and Data Analysis: Statistics for
    Censored Environmental Data, 1st edition. ed. Wiley-Interscience,
    Hoboken, N.J.
-   Helsel, D.R., 2011. Statistics for Censored Environmental Data Using
    Minitab© and R, 2nd edition. ed. Wiley-Interscience, Hoboken, N.J.

## Citing package <a name="cite"></a>

``` r
citation('NADA2')
```

    ## 
    ## To cite package 'NADA2' in publications use:
    ## 
    ##   Paul Julian and Dennis Helsel (2021). NADA2: Data Analysis for
    ##   Censored Environmental Data. R package version 1.0.0.
    ##   https://github.com/SwampThingPaul/NADA2
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Manual{,
    ##     title = {NADA2: Data Analysis for Censored Environmental Data},
    ##     author = {Paul Julian and Dennis Helsel},
    ##     year = {2021},
    ##     note = {R package version 1.0.0},
    ##     url = {https://github.com/SwampThingPaul/NADA2},
    ##   }

## Installation <a name="install"></a>

Development version can be installed from this repo.

``` r
install.packages("devtools");# if you do not have it installed on your PC
devtools::install_github("SwampThingPaul/NADA2")
```

------------------------------------------------------------------------
