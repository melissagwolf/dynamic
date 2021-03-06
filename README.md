
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dynamic

<!-- badges: start -->

<!-- badges: end -->

The goal of dynamic is to simulate fit index cutoffs for latent variable
models that are tailored to the user’s model statement, model type, and
sample size.

This is the counterpart of the Shiny Application,
[dynamicfit.app](https://dynamicfit.app/connect/). The Shiny app and the
R package will give you the same results. If you are comfortable with R,
consider using the package during high traffic times to reduce server
burden.

## Installation

This is the beta version of the package. Please submit bug reports and
issues on GitHub. You can install the beta version of dynamic from
[Github](https://github.com) with:

``` r
library(devtools)
devtools::install_github("melissagwolf/dynamic")
```

## Example

``` r
library(dynamic)

#Lavaan object example (manual=FALSE)
dat <- lavaan::HolzingerSwineford1939
lavmod <- "F1 =~ x1 + x2 + x3
          F2 =~ x4 + x5 + x6
          F3 =~ x7 + x8 + x9"
fit <- lavaan::cfa(lavmod,dat)
cfaHB(fit)

#Manual entry example (manual=TRUE)
manmod <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5 + .646*Y6"
cfaOne(manmod,500,manual=TRUE)
```
