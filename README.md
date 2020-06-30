
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dynamic

<!-- badges: start -->

<!-- badges: end -->

The goal of dynamic is to simulate model fit index cutoffs for latent
variable models that are tailored to the user’s model statement, model
type, and sample size. This is the counterpart of the Shiny Application,
[dynamicfit.app](https://dynamicfit.app/cfa). The Shiny app and the R
package will give you the same results. If you are comfortable with R,
consider using the package during high traffic times to reduce server
burden.

## Installation

This is the beta version of the package. Please submit bug reports and
issues on GitHub. You can install the beta version of dynamic from
[Github](https://github.com) with:

You can install the released version of dynamic from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("dynamic")
```

## Example

``` r
library(dynamic)
mod <- "F1 =~ .602*Y1 + .805*Y2 + .516*Y3 + .857*Y4
        F2 =~ .413*Y5 + -.631*Y6
        F1 ~~ .443*F2
        Y4 ~~ .301*Y5"
cfaFit(mod,500)
```

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
