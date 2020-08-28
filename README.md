
# rdirdirgamma

<!-- badges: start -->
<!-- badges: end -->

This package contains optimized code to sample from Dirichlet distribution, a specific hierarchical model (Dirichlet-Dirichlet-Gamma), and to perform ABC.  
Sampling is performed using Rcpp / RcppGSL.  
It has been written as part of my PhD thesis.

## Disclaimer

Code is provided as an exercise, do not expect it to be bug-free.  
Feel free to adapt, modify and reuse the code!

## Installation

You can install the development package from the GitHub repository:

``` r
remotes::install_github("lgaborini/rdirdirgamma")
```

### Requirements

This package requires GSL libraries.  
On Linux/Mac it should work without problems.

On Windows, one needs to install them separately:

1. download the `local323.zip` package from [here](http://www.stats.ox.ac.uk/pub/Rtools/libs.html) and extract it somewhere.
2. create the environmental variable `LIB_GSL` and point to the extracted path (e.g. `C:/local323_LIB_GSL`).  
You can do it using the `.Renviron` file, either project-wide or user-wide (`%HOME%/.Renviron`).
3. be sure that the file `src/Makevars.win` reads the GSL settings:

```sh
PKG_CPPFLAGS=-I$(LIB_GSL)/include -I../inst/include
PKG_LIBS=-L$(LIB_GSL)/lib/x64 -lgsl -lgslcblas
```

## Contents

### Dirichlet sampling

Sample from a Dirichlet distribution:

```r

library(rdirdirgamma)

alpha <- c(1, 1, 1)

x <- rdirichlet_cpp(alpha = alpha)
x
```

Sample from a Dirichlet distribution using the stick-breaking definition:

```r

library(rdirdirgamma)

n_samples <- 10
alpha <- c(1, 1, 1)

x <- rdirichlet_beta_cpp(n = n_samples, alpha = alpha)
x
```
