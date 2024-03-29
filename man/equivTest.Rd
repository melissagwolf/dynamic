% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equivTest.R
\name{equivTest}
\alias{equivTest}
\alias{print.equivTest}
\title{Equivalence testing with adjusted fit indexes for structural equation modeling}
\usage{
equivTest(
  n,
  T_ml = NULL,
  df = NULL,
  T_mli = NULL,
  p = NULL,
  manual = FALSE,
  plot = FALSE
)

\method{print}{equivTest}(x, ...)
}
\arguments{
\item{n}{This can either be a \code{\link{lavaan}} object, OR your sample size.}

\item{T_ml}{If you entered a \code{\link{lavaan}} object for n, leave this blank. Otherwise,
enter your model chi-square.}

\item{df}{If you entered a \code{\link{lavaan}} object for n, leave this blank. Otherwise,
enter your model degrees of freedom.}

\item{T_mli}{If you entered a \code{\link{lavaan}} object for n, leave this blank. Otherwise,
enter your baseline chi-square.}

\item{p}{If you entered a \code{\link{lavaan}} object for n, leave this blank. Otherwise,
enter the number of observed variables in your model.}

\item{manual}{If you entered a \code{\link{lavaan}} object, keep this set to FALSE.
If you manually entered each argument, set this to TRUE.}

\item{plot}{Displays a simple plot that compares your T-size RMSEA and T-Size CFI to the adjusted
bins.}

\item{x}{equivTest object}

\item{...}{other print parameters}
}
\value{
T-size RMSEA and T-Size CFI, along with adjusted bins for each index
}
\description{
This function generates adjusted fit index cutoffs using equivalence testing,
introduced by Yuan, Chan, Marcoulides, & Bentler (2016).
The default argument is a singular argument: a \code{\link{lavaan}} object.
The function can also accommodate manual entry of the sample size (n), model chi-square (T_ml),
degrees of freedom (df), baseline chi-square (T_mli), and number of observed variables (p).

The app-based version of this function can be found at \href{https://dynamicfit.app/}{dynamicfit.app}.
}
\examples{
#Lavaan object example (manual=FALSE)
dat <- lavaan::HolzingerSwineford1939
lavmod <- "F1 =~ x1 + x2 + x3
           F2 =~ x4 + x5 + x6
           F3 =~ x7 + x8 + x9"
fit <- lavaan::cfa(lavmod,dat)
equivTest(fit)

#' #Manual entry example (manual=TRUE)
n <- 301
T_ml <- 85.306
df <- 24
T_mli <- 918.852
p <- 9
equivTest(n,T_ml,df,T_mli,p,manual=TRUE)

}
\author{
Melissa G Wolf & Daniel McNeish

Maintainer: Melissa G Wolf <melissagordon@ucsb.edu>
}
