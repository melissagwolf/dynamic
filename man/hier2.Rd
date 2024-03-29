% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hier2.R
\name{hier2}
\alias{hier2}
\alias{print.hier2}
\title{Dynamic fit index (DFI) cutoffs for structural misspecification in hierarchical factor models}
\usage{
hier2(
  model,
  n = NULL,
  estimator = "ML",
  plot = FALSE,
  manual = FALSE,
  reps = 500
)

\method{print}{hier2}(x, ...)
}
\arguments{
\item{model}{This can either be a \code{\link{lavaan}} object from the \code{\link{cfa}} function,
OR a model statement written in \code{\link{lavaan}} \code{\link{model.syntax}} with standardized loadings.}

\item{n}{If you entered a \code{\link{lavaan}} object for model, leave this blank.
Otherwise, enter your sample size (numeric).}

\item{estimator}{Which estimator to use within the simulations (enter in quotes). The default is maximum likelihood.}

\item{plot}{Displays distributions of fit indices for each level of misspecification.}

\item{manual}{If you entered a \code{\link{lavaan}} object, keep this set to FALSE.
If you manually entered standardized loadings and sample size, set this to TRUE.}

\item{reps}{The number of replications used in your simulation. This is set to 500 by default in both the
R package and the corresponding Shiny App.}

\item{x}{hier2 object}

\item{...}{other print parameters}
}
\value{
Dynamic fit index (DFI) cutoffs for SRMR, RMSEA, and CFI.
}
\description{
This function generates DFI cutoffs for structural misspecification in hierarchical factor models using, by default, ML estimation.
The default argument is a singular argument: a \code{\link{lavaan}} object from the \code{\link{cfa}} function.
The function can also accommodate manual entry of the model statement and sample size.

The app-based version of this function can be found at \href{https://dynamicfit.app/}{dynamicfit.app}.
}
\examples{
#Manual entry example for a sample size of 2200 (manual=TRUE), from Reynolds & Keith (2017)

manmod <- "G =~ .51*F1 + .83*F2 + .97*F3 + .83*F4 + .87*F5 + .55*Y7
           F1 =~ .41*Y1 + .81*Y2 + .71*Y3
           F2 =~ .79*Y4 + .64*Y5 + .81*Y6 + .22*Y7
           F3 =~ .53*Y8 + .68*Y9 + .66*Y10
           F4 =~ .79*Y11 + .76*Y12
           F5 =~ .82*Y13 + .71*Y14 + .85*Y15 + .81*Y16
           F3 ~~ .77*F4"
\donttest{hier2(model=manmod,n=2200,manual=TRUE)}

}
\author{
Daniel McNeish, Melissa G Wolf, & Patrick D Manapat

Maintainer: Daniel McNeish <dmcneish@asu.edu>
}
