#' @title Reliability Representativeness of Coefficient Alpha or Omega
#'
#' @description This function evaluates how well a reliability summary index like alpha or omega
#'  represents the conditional reliability of a distribution of composite scores. It compares the conditional
#'  reliability function to a summary index and outputs a representativeness plot, a table of representativeness indices,and
#'  the full conditional reliability table for each possible sum score.
#'
#' @param data The original data to which the model was applied.
#' @param items Column names of the items on the scale being evaluated (entered as strings). If omitted, all variables in the data will be used.
#' @param rel Reliability coefficient to analyze. Options are "alpha" (the default) or "omega".
#' @param missing The missing data indicator in the data. Not needed in R, only present to simply use of this function in a Shiny application.
#' @param method how the test interval is created. Options are "CI" (the default), "width", or "raw".
#' "CI" uses a 95% Bayesian highest posterior density credible interval.
#' "width" builds an interval using a predetermined relative distance from the reliability coefficient (e.g., .05 from alpha).
#' "raw" builds an interval using a predetermined raw values (e.g., .70 to .90)
#' @param width Only required if method="width". Specifies a predetermined relative distance from the coefficient to each bound of the interval.
#' The total width of the interval will be twice this value
#' (e.g., if .05 is entered, the total interval width is .10 because it will span .05 above the coefficient and .05 below the coefficient)
#' @param raw.low Only required if method="raw". Manually specifies the lower bound of the test interval. Must be between 0 and 1.
#' @param raw.high Only required if method="raw". Manually specifies the upper bound of the test interval. Must be between 0 and 1.
#'
#' @import mirt ggplot2 Bayesrel stringr tidyr dplyr
#'
#' @author Daniel McNeish & Denis Dumas
#'
#' Maintainer: Daniel McNeish <dmcneish@asu.edu>
#'
#' @rdname RelRep
#'
#' @return Conditional reliability table and Reliability Representativeness plot and table.
#' @export
#'
#' @examples
#'
#' # "Example" dataset has 12 items on a 1-5 Likert scale
#'
#' #Example using the first 8 items in "Example" for testing coefficient alpha with a 95% Bayes credible interval
#' ex1<-RelRep(data=Example, items=c(names(Example[,1:8]), rel="alpha", method="CI"))
#'
#' #Example using odd items in "Example" to build a test interval that is .05 above and below coefficient omega
#' ex2<-RelRep(data=Example, items=c("X1","X3","X5","X7","X9","X11"), rel="omega", method="width", width=.05)
#'
#'#Example using even items in "Example" to specify a test interval for coefficient alpha between 0.70 and 0.80
#' ex3<-RelRep(data=Example, items=c("X2","X4","X6","X8","X10","X12"), rel="alpha", method="raw", raw.low=.70, raw.high=.80)

RelRep <- function(data, items=c(names(data)), rel="alpha", missing="NA", method="CI", width=NULL, raw.low=NULL, raw.high=NULL)
{
  if ( !(method%in% c("CI", "width", "raw"))){
    stop("Error:  the method input must be one of the following options: 'CI','width', or 'raw'")
    }

  if (method=="width" & is.null(width)){
    stop("Error: method='width' requires a user-defined value in the 'width=' option")
  }

  if (method=="raw" & (is.null(raw.low)|is.null(raw.high))){
    stop("Error: method='raw' requires a user-defined values both the 'raw.low=' and 'raw.high' options")
  }

  ##recode missing data to NA for R
  ##only used in Shiny app
  data[data == missing] <- NA

  #subset the uploaded data to only include selected items
  dat<-data[,items]

  ## delete any rows that are all missing to avoid errors
  dat<-dat[rowSums(is.na(dat)) != ncol(dat), ]

  #fit IRT model for which sum scores are a sufficient statistic
  model<-mirt::mirt(data=dat, itemtype="Rasch", verbose=F)

  ############
  #functions for conditional reliability
  #lifted from ggmirt package (Phillip Masurp, https://github.com/masurp/ggmirt)
  ############
  computeItemtrace <- function(pars, Theta, itemloc, offterm = matrix(0L, 1L, length(itemloc)-1L),
                               CUSTOM.IND, pis = NULL){
    if(is.null(pis)){
      itemtrace <- .Call('computeItemTrace', pars, Theta, itemloc, offterm)
      if(length(CUSTOM.IND)){
        for(i in CUSTOM.IND)
          itemtrace[,itemloc[i]:(itemloc[i+1L] - 1L)] <- ProbTrace(pars[[i]], Theta=Theta)
      }
    } else {
      tmp_itemtrace <- vector('list', length(pis))
      for(g in seq_len(length(pis))){
        tmp_itemtrace[[g]] <- .Call('computeItemTrace', pars[[g]]@ParObjects$pars, Theta, itemloc, offterm)
        if(length(CUSTOM.IND)){
          for(i in CUSTOM.IND)
            tmp_itemtrace[[g]][,itemloc[i]:(itemloc[i+1L] - 1L)] <- ProbTrace(pars[[g]]@ParObjects$pars[[i]], Theta=Theta)
        }
      }
      itemtrace <- do.call(rbind, tmp_itemtrace)
    }
    return(itemtrace)
  }

  ExtractGroupPars <- function(x){
    if(x@itemclass < 0L) return(list(gmeans=0, gcov=matrix(1)))
    nfact <- x@nfact
    gmeans <- x@par[seq_len(nfact)]
    phi_matches <- grepl("PHI", x@parnames)
    if (x@dentype == "Davidian") {
      phi <- x@par[phi_matches]
      tmp <- x@par[-c(seq_len(nfact), which(phi_matches))]
      gcov <- matrix(0, nfact, nfact)
      gcov[lower.tri(gcov, diag=TRUE)] <- tmp
      gcov <- makeSymMat(gcov)
      return(list(gmeans=gmeans, gcov=gcov, phi=phi))
    } else {
      par <- x@par
      if(x@dentype == "mixture") par <- par[-length(par)] # drop pi
      tmp <- par[-seq_len(nfact)]
      gcov <- matrix(0, nfact, nfact)
      gcov[lower.tri(gcov, diag=TRUE)] <- tmp
      gcov <- makeSymMat(gcov)
      return(list(gmeans=gmeans, gcov=gcov))
    }
  }

  makeSymMat <- function(mat){
    if(ncol(mat) > 1L){
      mat[is.na(mat)] <- 0
      mat <- mat + t(mat) - diag(diag(mat))
    }
    mat
  }

  #compute number of sum scores and number of people with each scores
  f<-mirt::fscores(model,method="EAPsum", full.scores=F,verbose=F, na.rm=TRUE)

  ## functions for computing conditional reliabilty (partially lifted for ggmirt package)
  nfact <- model@Model$nfact
  J <- model@Data$nitems
  theta <- f$F1
  ThetaFull <- Theta <- thetaComb(theta, nfact)
  info <- testinfo(model, ThetaFull)
  itemtrace <- computeItemtrace(model@ParObjects$pars, ThetaFull, model@Model$itemloc,
                                CUSTOM.IND=model@Internals$CUSTOM.IND)
  mins <- model@Data$mins
  maxs <- extract.mirt(model, 'K') + mins - 1
  gp <- ExtractGroupPars(model@ParObjects$pars[[J+1]])
  score <- c()
  for(i in 1:J)
    score <- c(score, (0:(model@Data$K[i]-1) + mins[i]) * (i %in% c(1:J)))
  score <- matrix(score, nrow(itemtrace), ncol(itemtrace), byrow = TRUE)
  plt <- data.frame(cbind(info,score=rowSums(score*itemtrace),Theta=Theta))
  colnames(plt) <- c("info", "score", "Theta")
  plt$SE <- 1 / sqrt(plt$info)
  plt$rxx <- plt$info / (plt$info + 1/gp$gcov[1L,1L])

  #possible sum scores in data
  plt$sum<-f$Sum.Scores
  #number of people with each score
  plt$observed<-f$observed

  if(method=="CI"){
    #compute alpha and 95% CI
    ci<-Bayesrel::strel(data=dat, estimates=rel, n.burnin=5000,n.iter=10000, missing="listwise")
    low<-as.numeric(ci$Bayes$cred$low)
    upp<-as.numeric(ci$Bayes$cred$up)
    est<-as.numeric(ci$Bayes$est)
  }

  if(method=="width"){
    ci<-Bayesrel::strel(data=dat, estimates=rel, n.burnin=5000,n.iter=10000, missing="listwise")
    est<-as.numeric(ci$Bayes$est)

    low<-est-width
    if(low<0){
      low<-0
        warning(
        "Warning: The 'width =' value implies an interval lower bound below 0.
         The lower bound has been replaced with 0.
         Check that your width was accurately entered", immediate.=T)
      }

    upp<-est+width

    if (upp>1) {
      upp<-1

      warning(
      "Warning: The 'width =' value implies an interval upper bound above 1.
       The upper bound has been replaced with 1.
       Check that your width was accurately entered", immediate.=T)
    }

  }

  if(method=="raw"){

    low<-raw.low
    upp<-raw.high

    if(upp>1) {
      stop(
        "Error: The 'raw.high =' option exceeds 1.
             Check that your width was accurately entered")
    }

    if (low<0) {
      stop(
        "Error: The 'low =' option is below 0.
             Check that your width was accurately entered")
    }

    if ((raw.low>raw.high | raw.low==raw.high)) {
      stop(
        "Error: The interval lower bound exceeds the interval upper bound.
                 Check that your width was accurately entered")
    }

    ci<-Bayesrel::strel(data=dat, estimates=rel, n.burnin=5000,n.iter=10000, missing="listwise")
    est<-as.numeric(ci$Bayes$est)
  }

  #percentage within CI, above CI, and below CI
  plt$within<-ifelse(low < plt$rxx & plt$rxx < upp, 1,0)
  plt$above<-ifelse(upp < plt$rxx, 1,0)
  plt$below<-ifelse(plt$rxx < low, 1,0)
  people<-round(100*apply(plt[,c("within","above","below")], 2, weighted.mean, w=plt$observed),2)

  #used for plot legend
  alpha <- data.frame(yintercept=est, Lines=str_to_title(rel))
  if(method=="CI"){
    CI <- data.frame(yintercept=low, Lines=paste(str_to_title(rel),'95% CI'))
  }

  if(method!="CI"){
    CI <- data.frame(yintercept=low, Lines='Test Interval')
  }


  #plot
  #conditional reliability, weighted by number of responses
  #with alpha and 95% CI
  #x axis is sum score, not theta score

  p<-ggplot(plt, aes(x = sum, y = rxx)) +
    geom_line(aes(color = (observed+1)),lwd=1.5) +
    scale_colour_gradient("People at Score",high ="red",low="blue")+
    ylim(0, 1) +
    theme_minimal(base_size=16)+
    geom_hline(yintercept=upp, linetype="dotted", lwd=.5)+
    geom_hline(aes(yintercept=yintercept,  linetype=Lines), CI, lwd=.5)+
    geom_hline(aes(yintercept=yintercept,  linetype=Lines), alpha, lwd=.5)+
    xlab("Sum Score") +
    ylab(paste("Relability"))


  #print alpha and percentages in console

  base::cat(paste0("Coefficient ", str_to_title(rel),":                  ", format(round(est,2),nsmall=2)))
  base::cat("\n")

  base::cat(paste0("Coefficient ", str_to_title(rel), " 95% CI:          ", "[",round(low,2),",",round(upp,2),"]"))
  base::cat("\n")
  base::cat("\n")

  base::cat(paste0("People above ", str_to_title(rel)," interval:    "))
  base::cat(paste0(unname(people[2]),"%"))
  base::cat("\n")

  base::cat(paste0("People within ", str_to_title(rel)," interval   :    "))
  base::cat(paste0(unname(people[1]),"%"))
  base::cat("\n")

  base::cat(paste0("People below  ", str_to_title(rel), " interval:    "))
  base::cat(paste0(unname(people[3]),"%"))


  ##Save a table for Shiny
  if(method=="CI"){
    y<-c(paste0("Coefficient ",str_to_title(rel),":"),paste0("Coefficient ", str_to_title(rel), " 95% CI:"), " ",
         paste0("People above ", str_to_title(rel)," interval:"),
         paste0("People within ", str_to_title(rel)," interval:"),
         paste0("People below ", str_to_title(rel)," interval :"),
         format(round(est,2),nsmall=2),
         paste0("[",format(round(low,2),nsmall=2),",",format(round(upp,2),nsmall=2),"]"),
         " ",
         paste0(unname(people[2]),"%"),
         paste0(unname(people[1]),"%"),
         paste0(unname(people[3]),"%"))
  }

  if(method!="CI"){
    y<-c(paste0("Coefficient ",str_to_title(rel),":"),paste0("Specified test interval:"), " ",
         paste0("People above ", str_to_title(rel)," interval:"),
         paste0("People within ", str_to_title(rel)," interval:"),
         paste0("People below ", str_to_title(rel)," interval :"),
         format(round(est,2),nsmall=2),
         paste0("[",format(round(low,2),nsmall=2),",",format(round(upp,2),nsmall=2),"]"),
         " ",
         paste0(unname(people[2]),"%"),
         paste0(unname(people[1]),"%"),
         paste0(unname(people[3]),"%"))
  }


  x<-matrix(y, nrow=6, ncol=2)

  x1<-data.frame(x)

  d<-plt[,c(6,5,7)]
  d[,2]<-round(d[,2],3)
  names(d)<-c("Sum Score","Conditional Reliability","Count")

  z<-list()

  #save summary table
  z$stat<-x1

  #save plot
  z$plot<-p

  #save conditional reliability table
  z$table<-d

  #print the plot by default
  print(z$p)

  #return the plot with the function
  return(z)

}
