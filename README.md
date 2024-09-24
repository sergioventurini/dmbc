# dmbc

<!-- badges: start -->
[![CRAN
status](https://www.r-pkg.org/badges/version/dmbc)](https://cran.r-project.org/package=dmbc)
[![](http://cranlogs.r-pkg.org/badges/grand-total/dmbc?color=blue)](https://cran.r-project.org/package=dmbc)
[![R-CMD-check](https://github.com/sergioventurini/dmbc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sergioventurini/dmbc/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

###### Current release: 1.0.3
###### R version required: at least 3.6.0
`R` package for Bayesian model-based clustering of several dissimilarity
matrices.

## Installation

Since the package requires some code to be compiled, you need a working C++
compiler. To get it:

- On Windows, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
- On Mac, install Xcode from the app store.
- On Linux, `sudo apt-get install r-base-dev` or similar.

Then, the easiest way to get the package is to install it from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("sergioventurini/dmbc")
```

See the help pages of the `dmbc()` and `dmbc_IC()` functions for examples
or run any of the demos available in the package by executing the code
`demo(_demo_name_, package = "dmbc")`.
