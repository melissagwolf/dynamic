#'@title DFI cutoffs for a Test of Exact Fit
#'
#'@description This function generates DFI cutoffs by treating the data generating model as the true model (using ML estimation).
#' The default argument is a singular argument: a \code{\link{lavaan}} object from the \code{\link{cfa}} function.
#' The function can also accommodate manual entry of the model statement and sample size.
#'
#' @param model This can either be a \code{\link{lavaan}} object from the \code{\link{cfa}} function,
#' OR a model statement written in \code{\link{lavaan}} \code{\link{model.syntax}} with standardized loadings.
#' @param n If you entered a \code{\link{lavaan}} object for model, leave this blank.
#' Otherwise, enter your sample size (numeric).
#' @param reps The number of replications used in your simulation (default is 1000).
#' @param plot Displays distributions of fit indices for each fit index.
#' @param manual If you entered a \code{\link{lavaan}} object, keep this set to FALSE.
#' If you manually entered standardized loadings and sample size, set this to TRUE.
#'
#' @import dplyr lavaan simstandard
#' @importFrom purrr map map_dfr
#' @importFrom tidyr nest
#' @importFrom graphics hist abline
#'
#' @author Melissa G Wolf & Daniel McNeish
#'
#' Maintainer: Melissa G Wolf <melissagordon@ucsb.edu>
#'
#' @rdname exactFit
#'
#' @return Dynamic fit index (DFI) cutoffs for Chi-Square, SRMR, RMSEA, and CFI.
#' @export
#'
#' @examples
#' #Lavaan object example (manual=FALSE)
#' dat <- lavaan::HolzingerSwineford1939
#' lavmod <- "F1 =~ x1 + x2 + x3
#'            F2 =~ x4 + x5 + x6
#'            F3 =~ x7 + x8 + x9"
#' fit <- lavaan::cfa(lavmod,dat)
#' exactFit(fit)
#'
#' #Manual entry example (manual=TRUE)
#' manmod <- "F1 =~ .602*Y1 + .805*Y2 + .516*Y3 + .857*Y4
#'            F2 =~ .413*Y5 + -.631*Y6
#'            F1 ~~ .443*F2
#'            Y4 ~~ .301*Y5"
#' exactFit(manmod,500,manual=TRUE)

exactFit <- function(model,n,reps=1000,plot=FALSE,manual=FALSE){

  #If manual, expect manual (a la Shiny app)
  if(manual){
    model9 <- model
    n <- n
  }else{
    #Use this to rewrite error message for when someone forgot to use manual=TRUE
    #But entered in model statement and sample size
    #https://community.rstudio.com/t/create-custom-error-messages/39058/4
    #This is hacky but works, although traceback might confuse people
    tryCatch(cfa_n(model),
             error=function(err){
               if (grepl("trying to get slot", err)) {
                 stop("dynamic Error: Did you forget to use manual=TRUE?")
               }
             })

    #Error for when someone enters an object that doesn't exist, or a non-lavaan object
    tryCatch(cfa_n(model),
             error=function(err2){
               if (grepl("Error in base::unlist", err2)){
                 stop("dynamic Error: Did you enter a lavaan object? Confirm that it is a lavaan object using class(). If you do not have a lavaan object, enter the arguments manually and select manual=TRUE.")
               }
             })

    #Use these functions to convert to manual (input is a lavaan object)
    #Probably what we should expect for people using R
    #need 'n' first because otherwise model will overwrite
    n <- cfa_n(model)
    model9 <- cfa_lavmod(model)

  }

  if (unstandardized(model9)>0){
    stop("dynamic Error: One of your loadings or correlations has an absolute value of 1 or above (an impossible value). Please use standardized loadings. If all of your loadings are under 1, try looking for a missing decimal somewhere in your model statement.")
  }

  #Create list to store outputs (table and plot)
  res <- list()

  #Output fit indices if someone used manual=F
  #Will ignore in print statement if manual=T
  #Exclamation point is how we indicate if manual = T (because default is F)

  if(!manual){
    fitted <- round(lavaan::fitmeasures(model,c("chisq","df","pvalue","srmr","rmsea","cfi")),3)
    fitted_m <- as.matrix(fitted)
    fitted_t <- t(fitted_m)
    fitted_t <- as.data.frame(fitted_t)
    colnames(fitted_t) <- c("Chi-Square"," df","p-value","  SRMR","  RMSEA","   CFI")
    rownames(fitted_t) <- c("")
    res$fit <- fitted_t
  }

  #Create object type (S4)
  #setClass("res",slots=list(dat="data.frame",cutoffs="data.frame",plots="list"))
  #https://www.datacamp.com/community/tutorials/r-objects-and-classes

  #Run simulation
  exact_dat <- exact_fit_dat(model9,n,reps)

  #Save the data
  res$data <- exact_dat

  #Extract cutoff values
  exact_vals <- exact_dat %>%
    dplyr::summarise(chisq=round(quantile(chisq,c(.99,.95)),3),
                     df=mean(df),
                     srmr=round(quantile(srmr,c(.99,.95)),3),
                     rmsea=round(quantile(rmsea,c(.99,.95)),3),
                     cfi=round(quantile(cfi,c(.01,.05)),3))

  #row names for tibbles is deprecated - might need to convert to df in the future
  exact_vals <- as.data.frame(exact_vals)
  colnames(exact_vals) <- c("Chi-Square"," df","  SRMR","  RMSEA","   CFI")
  rownames(exact_vals) <- c("99th:","95th:")

  #Put into output list
  res$cutoffs <- exact_vals

  #plots=T

  if(plot){
    #Create basic histograms
    fig_chi <- graphics::hist(exact_dat$chisq,
                    main="Chi-Square Histogram",
                    xlab="Chi-Square",
                    breaks="FD")
    graphics::abline(v=exact_vals$`Chi-Square`[1],col="red")

    fig_srmr <- graphics::hist(exact_dat$srmr,
                     main="SRMR Histogram",
                     xlab="SRMR",
                     breaks="FD")
    graphics::abline(v=exact_vals$`  SRMR`[1],col="blue")

    fig_rmsea <- graphics::hist(exact_dat$rmsea,
                      main="RMSEA Histogram",
                      xlab="RMSEA",
                      breaks="FD")
    graphics::abline(v=exact_vals$`  RMSEA`[1],col="purple")

    fig_cfi <- graphics::hist(exact_dat$cfi,
                    main="CFI Histogram",
                    xlab="CFI",
                    breaks="FD")
    graphics::abline(v=exact_vals$`   CFI`[1],col="green")

    plots <- list(fig_chi,fig_srmr,fig_rmsea,fig_cfi)

    res$plots <- plots
  }

  #Create object (necessary for subsequent print statement)
  class(res) <- 'exactFit'

  return(res)

}

#' @method print exactFit
#' @param x exactFit object
#' @param ... other print parameters
#' @rdname exactFit
#' @export

print.exactFit <- function(x,...){

  #Automatically return fit cutoffs
  base::cat("DFI cutoffs: \n")
  base::print(x$cutoffs)

  #Only print fit indices from lavaan object if someone submits a lavaan object
  if(!is.null(x$fit)){
    base::cat("\n")

    base::cat("Empirical fit indices: \n")
    base::print(x$fit)
  }

  #Currently automatically returning plots when T, which is annoying
  #Prefer them hidden in the list unless called

  #Hides this function
  base::invisible()
}
