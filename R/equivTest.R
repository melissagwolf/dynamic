#' @title Equivalence testing with adjusted fit indexes for structural equation modeling
#'
#' @description This function generates adjusted fit index cutoffs using equivalence testing,
#' introduced by Yuan, Chan, Marcoulides, & Bentler (2016).
#' The default argument is a singular argument: a \code{\link{lavaan}} object.
#' The function can also accommodate manual entry of the sample size (n), model chi-square (T_ml),
#' degrees of freedom (df), baseline chi-square (T_mli), and number of observed variables (p).
#'
#' The app-based version of this function can be found at \href{https://dynamicfit.app/}{dynamicfit.app}.
#'
#' @param n This can either be a \code{\link{lavaan}} object, OR your sample size.
#' @param T_ml If you entered a \code{\link{lavaan}} object for n, leave this blank. Otherwise,
#' enter your model chi-square.
#' @param df If you entered a \code{\link{lavaan}} object for n, leave this blank. Otherwise,
#' enter your model degrees of freedom.
#' @param T_mli If you entered a \code{\link{lavaan}} object for n, leave this blank. Otherwise,
#' enter your baseline chi-square.
#' @param p If you entered a \code{\link{lavaan}} object for n, leave this blank. Otherwise,
#' enter the number of observed variables in your model.
#' @param manual If you entered a \code{\link{lavaan}} object, keep this set to FALSE.
#' If you manually entered each argument, set this to TRUE.
#' @param plot Displays a simple plot that compares your T-size RMSEA and T-Size CFI to the adjusted
#' bins.
#'
#' @import stringr ggplot2
#'
#' @author Melissa G Wolf & Daniel McNeish
#'
#' Maintainer: Melissa G Wolf <melissagordon@ucsb.edu>
#'
#' @rdname equivTest
#'
#' @return T-size RMSEA and T-Size CFI, along with adjusted bins for each index
#' @export
#'
#' @examples
#' #Lavaan object example (manual=FALSE)
#' dat <- lavaan::HolzingerSwineford1939
#' lavmod <- "F1 =~ x1 + x2 + x3
#'            F2 =~ x4 + x5 + x6
#'            F3 =~ x7 + x8 + x9"
#' fit <- lavaan::cfa(lavmod,dat)
#' equivTest(fit)
equivTest <- function(n,T_ml=NULL,df=NULL,T_mli=NULL,p=NULL,manual=FALSE,plot=FALSE){

  #if manual, expect manual (a la shiny app)
  if(manual){
    p=p
    T_ml=T_ml
    df=df
    T_mli=T_mli
    n9=n
  }else{
    #Use this to rewrite error message for when someone forgot to use manual=TRUE
    #But entered in model statement and sample size
    #This is hacky but works, although traceback might confuse people
    #https://community.rstudio.com/t/create-custom-error-messages/39058/4
    tryCatch(equiv_T_ml(n),
             error=function(err){
               if (grepl("trying to get slot", err)) {
                 stop("dynamic Error: Did you forget to use manual=TRUE?")
               }
             })

    #Error for when someone enters an object that doesn't exist, or a non-lavaan object
    tryCatch(equiv_T_ml(n),
             error=function(err2){
               if (grepl("Error in base::unlist", err2)){
                 stop("dynamic Error: Did you enter a lavaan object? Confirm that it is a lavaan object using class(). If you do not have a lavaan object, enter the arguments manually and select manual=TRUE.")
               }
             })

    #Use these functions to convert to manual (input is a lavaan object)
    #Probably what we should expect for people using R
    #need 'n' last because otherwise model will overwrite
    #no longer need n last with n9
    T_ml <- equiv_T_ml(n)
    df <- equiv_df(n)
    T_mli <- equiv_T_mli(n)
    p <- equiv_p(n)
    n9 <- equiv_n(n)
  }

  #Error for when someone doesn't enter all 5 arguments when they select manual=T
  tryCatch(equiv_cutoffs(p,T_ml,df,T_mli,n9),
           error=function(err3){
             if (grepl("non-numeric argument to binary operator", err3)){
               stop("dynamic Error: Did you enter all 5 arguments and select manual=TRUE?")
             }
           })

  #Create list to store outputs (table and plot)
  res <- list()

  #Output fit indices if someone used manual=F
  #Will ignore in print statement if manual=T
  #Exclamation point is how we indicate if manual = T (because default is F)

  if(!manual){
    fitted <- round(lavaan::fitmeasures(n,c("chisq","df","pvalue","srmr","rmsea","cfi")),3)
    fitted_m <- as.matrix(fitted)
    fitted_t <- t(fitted_m)
    fitted_t <- as.data.frame(fitted_t)
    colnames(fitted_t) <- c("Chi-Square"," df","p-value","  SRMR","  RMSEA","   CFI")
    rownames(fitted_t) <- c("")
    res$fit <- fitted_t
  }

  #Get data in
  dat <- equiv_cutoffs(p,T_ml,df,T_mli,n9)

  #Remove 0's
  clean <- lapply(dat, function(x) stringr::str_replace_all(x,"0\\.","."))
  ul <- unlist(clean)
  cut <- matrix(ul,nrow=2,ncol=5)

  ## Extract T-Size and save to list
  rmsea <- cut[1,5]
  res$rmsea <- rmsea

  cfi <- cut[2,5]
  res$cfi <- cfi

  ## Label cutoffs
  excellent <- "Excellent: "
  close <- "Close: "
  fair <- "Fair: "
  mediocre <- "Mediocre: "
  poor <- "Poor: "

  ## Create RMSEA bins
  one_r <- paste(cut[1,1],"or below")
  two_r <- paste(cut[1,1],"to",cut[1,2])
  three_r <- paste(cut[1,2],"to",cut[1,3])
  four_r <- paste(cut[1,3],"to",cut[1,4])
  five_r <- paste(cut[1,4],"or above")

  #Combine
  eo_r <- paste(excellent,one_r,sep="")
  ct_r <- paste(close,two_r,sep="")
  ft_r <- paste(fair,three_r,sep="")
  mf_r <- paste(mediocre,four_r,sep="")
  pf_r <- paste(poor,five_r,sep="")

  #Save to list
  res$eo_r <- eo_r
  res$ct_r <- ct_r
  res$ft_r <- ft_r
  res$mf_r <- mf_r
  res$pf_r <- pf_r

  #Create CFI Bins
  one_c <- paste(cut[2,1],"or below")
  two_c <- paste(cut[2,1],"to",cut[2,2])
  three_c <- paste(cut[2,2],"to",cut[2,3])
  four_c <- paste(cut[2,3],"to",cut[2,4])
  five_c <- paste(cut[2,4],"or above")

  #Combine
  eo_c <- paste(poor,one_c,sep="")
  ct_c <- paste(mediocre,two_c,sep="")
  ft_c <- paste(fair,three_c,sep="")
  mf_c <- paste(close,four_c,sep="")
  pf_c <- paste(excellent,five_c,sep="")

  #Save to list
  res$eo_c <- eo_c
  res$ct_c <- ct_c
  res$ft_c <- ft_c
  res$mf_c <- mf_c
  res$pf_c <- pf_c

  if (plot){

    #RMSEA
    e <- max(dat[1,4],dat[1,5])

    x <- dat[1,1:4]

    m <- e+(dat[1,1]-0)

    ex <- mean(c(.00001,x[1]))
    cl <- mean(c(x[1],x[2]))
    fa <- mean(c(x[2],x[3]))
    me <- mean(c(x[3],x[4]))
    po <- mean(c(x[4],m))

    rmsea_plot <- ggplot(data.frame(x), aes(x=x, y=0)) +
      geom_point(alpha=0)  +
      annotate("segment",x=0,xend=m, y=0, yend=0, size=1,col="grey50") +
      annotate("segment",x=0,xend=0, y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=m,xend=m, y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[1],xend=x[1], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[2],xend=x[2], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[3],xend=x[3], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[4],xend=x[4], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=dat[1,5],xend=dat[1,5],y=-0.1,yend=.25, size=1, col="tomato4")+
      annotate("text",x=dat[1,5],y=.6,label=paste("T-size \n RMSEA \n",cut[1,5]),
               col="tomato4",size=3.5)+
      annotate("text",x=ex,y=-.5,label="Excellent",size=3.5)+
      annotate("text",x=cl,y=-.5,label="Close",size=3.5)+
      annotate("text",x=fa,y=-.5,label="Fair",size=3.5)+
      annotate("text",x=me,y=-.5,label="Mediocre",size=3.5)+
      annotate("text",x=po,y=-.5,label="Poor",size=3.5)+
      geom_text(aes(label = x),col="grey20", position=position_nudge(y=-.2),size=3.5) +
      scale_x_continuous(limits = c(0,m)) +
      scale_y_continuous(limits = c(-1,1)) +
      scale_color_manual(values = unname(colours)) +
      theme(panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())

    #CFI
    e <- min(dat[2,1],dat[2,5])

    x <- dat[2,1:4]

    m <- e-(1-dat[2,4])

    ex <- mean(c(1,x[4]))
    cl <- mean(c(x[4],x[3]))
    fa <- mean(c(x[3],x[2]))
    me <- mean(c(x[2],x[1]))
    po <- mean(c(x[1],m))

    cfi_plot <- ggplot(data.frame(x), aes(x=x, y=0)) +
      geom_point(alpha=0)  +
      annotate("segment",x=m,xend=1, y=0, yend=0, size=1,col="grey50") +
      annotate("segment",x=1,xend=1, y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=m,xend=m, y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[1],xend=x[1], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[2],xend=x[2], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[3],xend=x[3], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=x[4],xend=x[4], y=-0.1,yend=0.1, size=1,col="grey50") +
      annotate("segment",x=dat[2,5],xend=dat[2,5],y=-0.1,yend=.25, size=1, col="tomato4")+
      annotate("text",x=dat[2,5],y=.6,label=paste("T-size \n CFI \n",cut[2,5]),
               col="tomato4", size=3.5)+
      annotate("text",x=ex,y=-.5,label="Excellent",size=3.5)+
      annotate("text",x=cl,y=-.5,label="Close",size=3.5)+
      annotate("text",x=fa,y=-.5,label="Fair",size=3.5)+
      annotate("text",x=me,y=-.5,label="Mediocre",size=3.5)+
      annotate("text",x=po,y=-.5,label="Poor",size=3.5)+
      geom_text(aes(label = x),col="grey20", position=position_nudge(y=-.2),size=3.5) +
      scale_x_continuous(limits = c(m,1)) +
      scale_y_continuous(limits = c(-1,1)) +
      scale_color_manual(values = unname(colours)) +
      theme(panel.background = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank())

    res$R_Plot <- rmsea_plot
    res$C_Plot <- cfi_plot
  }

  class(res) <- 'equivTest'

  return(res)
}

#' @method print equivTest
#' @param x equivTest object
#' @param ... other print parameters
#' @rdname equivTest
#' @export

#Print suppression/organization statement for list
#Needs same name as class, not function name
#Need to add ... param or will get error message in CMD check

print.equivTest <- function(x,...){

  #Automatically return fit cutoffs
  cat("\n",
      "Your T-size RMSEA:",
      x$rmsea,
      "\n","\n",
      "Adjusted RMSEA cutoff values: \n",
      x$eo_r,"\n",
      x$ct_r,"\n",
      x$ft_r,"\n",
      x$mf_r,"\n",
      x$pf_r,
      "\n","\n",
      "Your T-size CFI:", x$cfi,
      "\n","\n",
      "Adjusted CFI cutoff values: \n",
      x$pf_c,"\n",
      x$mf_c,"\n",
      x$ft_c,"\n",
      x$ct_c,"\n",
      x$eo_c,"\n")

  #Only print fit indices from lavaan object if someone submits a lavaan object
  if(!is.null(x$fit)){
    base::cat("\n")

    base::cat("Empirical fit indices: \n")
    base::print(x$fit)
  }

  if(!is.null(x$R_Plot)){

    cat("\n",
        "The plots for each fit index are in the Plots tab")
    print(x$R_Plot)
    print(x$C_Plot)
  }

  #Hides this function
  base::invisible()
}
