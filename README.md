# dmbc

<!-- badges: start -->

[![Travis build status](https://travis-ci.org/sergioventurini/dmbc.svg?branch=master)](https://travis-ci.org/sergioventurini/dmbc)
[![CRAN
status](https://www.r-pkg.org/badges/version/dmbc)](https://cran.r-project.org/package=dmbc)
<!-- badges: end -->

## Overview

###### Current release: 0.2.3
###### R version required: at least 3.5.0
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

See the help pages of the `dmbc()` and `dmbc_IC()` functions for some examples
or have a look at the demos in the package.
