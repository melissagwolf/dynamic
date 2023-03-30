#' @title Dynamic fit index (DFI) cutoffs for continuous, non-normal one-factor CFA models with (possible) missing data
#'
#' @description This function generates DFI cutoffs for one-factor CFA models that treats items as continuous and non-normal with
#' possible missing data. This functions uses a modified Bollen-Stine bootstrap to accommodate non-normality and missingness rather
#' than simulating from a particular distribution. The default argument is a singular argument: a \code{\link{lavaan}} object from
#'  the \code{\link{cfa}} function. The function can also accommodate manual entry of the model statement and
#'  sample size (including threshold estimates). A primary difference in nnor DFI functions is that a dataset from which to bootstrap
#'  must also be provided in the 'data' argument.
#'
#' The app-based version of this function can be found at \href{https://dynamicfit.app/}{dynamicfit.app}.
#'
#' @param model This can either be a \code{\link{lavaan}} object from the \code{\link{cfa}} function,
#' OR a model statement written in \code{\link{lavaan}} \code{\link{model.syntax}} with standardized loadings and thresholds.
#' @param data An empirical dataset to which a modified Bollen-Stine bootstrap will be applied to create hypothetical misspecified data
#' @param n If you entered a \code{\link{lavaan}} object for model, leave this blank.
#' Otherwise, enter your sample size (numeric).
#' @param plot Displays distributions of fit indices for each level of misspecification.This also includes plots to visualize how close the
#' distributions of the hypothetical data come to the original data.
#' @param manual If you entered a \code{\link{lavaan}} object, keep this set to FALSE.
#' If you manually entered standardized loadings and sample size, set this to TRUE.
#' @param reps The number of replications used in your simulation. This is set to 500 by default in both the
#' R package and the corresponding Shiny App.
#' @param estimator Which estimator to use within the simulations (enter in quotes). The default is MLR
#'
#' @import dplyr lavaan simstandard ggplot2 stringr
#' @importFrom purrr map map_dfr map2
#' @importFrom tidyr unite extract
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#' @importFrom semTools bsBootMiss
#'
#' @author Daniel McNeish & Melissa G Wolf
#'
#' Maintainer: Daniel McNeish <dmcneish@asu.edu>
#'
#' @rdname nnorOne
#'
#' @return Dynamic fit index (DFI) cutoffs for SRMR, RMSEA, and CFI.
#' @export
#'
#' @examples
#' #Lavaan object example (manual=FALSE)
#' dat <- lavaan::HolzingerSwineford1939
#' lavmod <- "F1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9"
#' fit <- lavaan::cfa(lavmod,dat)
#' \donttest{cfaOne(fit)}
#'
#' #Manual entry example for a sample size of 6649 (manual=TRUE)
#' #thresholds go in model statement as (a)categorical item name, (b) vertical pipe, (c) estimate, (d)times t+threshold number
#' manmod <- "E =~ .76*bfi_e1 + .73*bfi_e2 + .59*bfi_e3 + .71*bfi_e4 + .84*bfi_e5 + .58*bfi_e6 + .71*bfi_e7 + .80*bfi_e8
#'
#'bfi_e1 |-1.69*t1
#'bfi_e1 |-1.06*t2
#'bfi_e1 |-0.53* t3
#'bfi_e1 |0.06* t4
#'bfi_e1 |0.75* t5
#'
#'bfi_e2 |-1.16*t1
#'bfi_e2 |-0.42*t2
#'bfi_e2 |0.28*t3
#'bfi_e2 |0.71*t4
#'bfi_e2 |1.34*t5
#'
#'bfi_e3 |-1.99* t1
#'bfi_e3 |-1.27* t2
#'bfi_e3 | 0.61*t3
#'bfi_e3 | 0.09*t4
#'bfi_e3 | 0.97*t5
#'
#'bfi_e4 |-2.05* t1
#'bfi_e4 | -1.36*t2
#'bfi_e4 | -0.74*t3
#'bfi_e4 | -.004*t4
#'bfi_e4 | 0.81*t5
#'
#'bfi_e5 | -1.12*t1
#'bfi_e5 | -0.42*t2
#'bfi_e5 | 0.16*t3
#'bfi_e5 | 0.60*t4
#'bfi_e5 | 1.21*t5
#'
#'bfi_e6 | -1.76*t1
#'bfi_e6 | -1.18*t2
#'bfi_e6 | -0.68*t3
#'bfi_e6 | -0.03*t4
#'bfi_e6 | 0.76*t5
#'
#'bfi_e7 | -1.04*t1
#'bfi_e7 | -0.31*t2
#'bfi_e7 | 0.46*t3
#'bfi_e7 | 0.78*t4
#'bfi_e7 | 1.38*t5
#'
#'bfi_e8 | -1.74*t1
#'bfi_e8 | -1.12*t2
#'bfi_e8 | -0.64*t3
#'bfi_e8 | -.002*t4
#'bfi_e8 | 0.77*t5"
#' \donttest{catOne(model=manmod,n=6649,manual=TRUE)}
#'

nnorOne <- function(model,data,n=NULL,plot=FALSE,manual=FALSE,estimator="MLR",reps=500){

  #If manual, expect manual (a la Shiny app)
  if(manual){
    n <- n
    model9 <- model
  }else{
    #Use this to rewrite error message for when someone forgot to use manual=TRUE
    #But entered in model statement and sample size
    #This is hacky but works, although traceback might confuse people
    #https://community.rstudio.com/t/create-custom-error-messages/39058/4
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

  if (number_factor(model9)>1){
    stop("dynamic Error: You entered a multi-factor model.  Use cfaHB instead.")
  }

  if (defre(model9,n)==0){
    stop("dynamic Error: It is impossible to add misspecifications to a just identified model.")
  }

  if ( nrow(one_num(model9)) < (number_factor(model9)-1)){
    stop("dynamic Error: There are not enough free items to produce all misspecification levels.")
  }

  #Create list to store outputs (table and plot)
  res <- list()

  #Output fit indices if someone used manual=F
  #Will ignore in print statement if manual=T
  #Exclamation point is how we indicate if manual = T (because default is F)

  if(!manual){
    if (model@Options$test=="satorra.bentler" |model@Options$test=="yuan.bentler.mplus" | model@Options$test=="yuan.bentler.mplus"){
      fitted <- round(lavaan::fitmeasures(model,c("chisq.scaled","df","pvalue.scaled","srmr","rmsea.robust","cfi.robust")),3)
    } else if (model@Options$test=="scaled.shifted" | model@Options$test=="mean.var.adusted"){
      fitted <- round(lavaan::fitmeasures(model,c("chisq.scaled","df","pvalue.scaled","srmr","rmsea.scaled","cfi.scaled")),3)
    } else if(model@Options$test=="standard" ){
      fitted <- round(lavaan::fitmeasures(model,c("chisq","df","pvalue","srmr","rmsea","cfi")),3)}
    fitted_m <- as.matrix(fitted)
    fitted_t <- t(fitted_m)
    fitted_t <- as.data.frame(fitted_t)
    colnames(fitted_t) <- c("Chi-Square"," df","p-value","  SRMR","  RMSEA","   CFI")
    rownames(fitted_t) <- c("")
    res$fit <- fitted_t
  }

  #Run simulation
  tryCatch(one_df_nnor(model9,data,n,estimator, reps),
           error=function(err3){
             if (grepl("is missing, with no default", err3)){
               stop("dynamic Error: Did you forget to include a dataset? nnor functions in dynamic use bootstrapping instead of simulation to ensure
                    that distributions and missing data patterns match the original data. The nnorOne function requires a dataset an input.")
             }
           })
  results <- one_df_nnor(model9,data,n,estimator, reps)

  #Save the data and make it exportable
  res$data <- fit_data(results)

  #For each list element (misspecification) compute the cutoffs
  misspec_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_M=stats::quantile(SRMR_M, c(seq(0.05,1,0.01))),
                                                      RMSEA_M=stats::quantile(RMSEA_M, c(seq(0.05,1,0.01))),
                                                      CFI_M=stats::quantile(CFI_M, c(seq(0.95,0,-0.01)))))

  #For the true model, compute the cutoffs (these will all be the same - just need in list form)
  true_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_T=stats::quantile(SRMR_T, c(.95)),
                                                   RMSEA_T=stats::quantile(RMSEA_T, c(.95)),
                                                   CFI_T=stats::quantile(CFI_T, c(.05))))

  fit<-list()
  S<-list()
  R<-list()
  C<-list()
  final<-list()
  for (i in 1:length(misspec_sum))
  {
    fit[[i]]<-cbind(misspec_sum[[i]], true_sum[[i]])
    fit[[i]]$Power<-seq(.95, 0.0, -.01)
    fit[[i]]$S<-ifelse(fit[[i]]$SRMR_M >= fit[[i]]$SRMR_T, 1, 0)
    fit[[i]]$R<-ifelse(fit[[i]]$RMSEA_M >= fit[[i]]$RMSEA_T, 1, 0)
    fit[[i]]$C<-ifelse(fit[[i]]$CFI_M <= fit[[i]]$CFI_T, 1, 0)
    S[[i]]<-subset(fit[[i]], subset=!duplicated(fit[[i]][('S')]), select=c("SRMR_M","Power","S")) %>% filter(S==1)
    R[[i]]<-subset(fit[[i]], subset=!duplicated(fit[[i]][('R')]), select=c("RMSEA_M","Power","R")) %>% filter(R==1)
    C[[i]]<-subset(fit[[i]], subset=!duplicated(fit[[i]][('C')]), select=c("CFI_M","Power","C"))  %>% filter(C==1)
    colnames(S[[i]])<-c("SRMR","PowerS")
    colnames(R[[i]])<-c("RMSEA","PowerR")
    colnames(C[[i]])<-c("CFI","PowerC")

    final[[i]]<-cbind(S[[i]],R[[i]],C[[i]])
    final[[i]]<-final[[i]][c("SRMR","PowerS","RMSEA","PowerR","CFI","PowerC")]
  }

  L0<-data.frame(cbind(true_sum[[1]]$SRMR_T,.95,true_sum[[1]]$RMSEA_T,0.95,true_sum[[1]]$CFI_T,0.95))%>%
    `colnames<-`(c("SRMR","PowerS","RMSEA","PowerR","CFI","PowerC"))

  if(length(misspec_sum)==3) {
    Fit<-round(rbind(L0,final[[1]],final[[2]],final[[3]]),3)
  }

  if(length(misspec_sum)==2) {
    Fit<-round(rbind(L0,final[[1]],final[[2]]),3)
  }

  if(length(misspec_sum)==1) {
    Fit<-round(rbind(L0,final[[1]]),3)
  }

  fit1<-unlist(Fit)%>% matrix(nrow=(length(misspec_sum)+1), ncol=6) %>%
    `colnames<-`(c("SRMR","PowerS","RMSEA","PowerR","CFI","PowerC"))

  PS<-paste(round(100*fit1[,2], 2), "%", sep="")
  PR<-paste(round(100*fit1[,4], 2), "%", sep="")
  PC<-paste(round(100*fit1[,6], 2), "%", sep="")

  for (j in 2:(length(misspec_sum)+1)) {
    if(fit1[j,2]<.50){fit1[j,1]<-"NONE"}
    if(fit1[j,4]<.50){fit1[j,3]<-"NONE"}
    if(fit1[j,6]<.50){fit1[j,5]<-"NONE"}
  }
  fit1[,2]<-PS
  fit1[,4]<-PR
  fit1[,6]<-PC

  pp<-c(rep("--",(length(misspec_sum)+1)))
  pp0<-c(rep("",(length(misspec_sum)+1)))

  SS<-noquote(matrix(rbind(fit1[,1],fit1[,2],pp0),ncol=1))
  RR<-noquote(matrix(rbind(fit1[,3],fit1[,4],pp0),ncol=1))
  CC<-noquote(matrix(rbind(fit1[,5],fit1[,6],pp0),ncol=1))

  Table<-noquote(cbind(SS,RR,CC) %>%
                   `colnames<-`(c("SRMR","RMSEA","CFI")))

  if(length(misspec_sum)==3) {
    rownames(Table)<-c("Level-0","Specificity", "","Level-1", "Sensitivity","", "Level-2", "Sensitivity","", "Level-3", "Sensitivity","")
  }

  if(length(misspec_sum)==2) {
    rownames(Table)<-c("Level-0","Specificity", "","Level-1", "Sensitivity","", "Level-2", "Sensitivity","")
  }

  if(length(misspec_sum)==1) {
    rownames(Table)<-c("Level-0","Specificity", "","Level-1", "Sensitivity","")
  }


  Table<-Table[1:(nrow(Table)-1),]

  #Put into list
  res$cutoffs <- Table

  #If user selects plot = T
  if(plot){

    #For each list element (misspecification) compute the cutoffs
    #misspec_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_M=stats::quantile(SRMR_M, c(.05,.1)),
    #                                                    RMSEA_M=stats::quantile(RMSEA_M, c(.05,.1)),
    #                                                   CFI_M=stats::quantile(CFI_M, c(.95,.9))))

    #For the true model, compute the cutoffs (these will all be the same - just need in list form)
    # true_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_T=stats::quantile(SRMR_T, c(.95,.9)),
    #                                                 RMSEA_T=stats::quantile(RMSEA_T, c(.95,.9)),
    #                                                 CFI_T=stats::quantile(CFI_T, c(.05,.1))))

    #Select just those variables and rename columns to be the same
    Misspec_dat <- purrr::map(results,~dplyr::select(.,SRMR_M:Type_M) %>%
                                `colnames<-`(c("SRMR","RMSEA","CFI","Model")))

    #Select just those variables and rename columns to be the same
    True_dat <- purrr::map(results,~dplyr::select(.,SRMR_T:Type_T) %>%
                             `colnames<-`(c("SRMR","RMSEA","CFI","Model")))

    #For each element in the list, bind the misspecified cutoffs to the true cutoffs
    #rbind doesn't work well with lists (needs do.call statement)
    plot <- base::lapply(base::seq(base::length(Misspec_dat)),function(x) dplyr::bind_rows(Misspec_dat[x],True_dat[x]))

    #Plot SRMR. Need map2 and data=.x (can't remember why).
    SRMR_plot <- purrr::map2(plot,final,~ggplot(data=.x,aes(x=SRMR,fill=Model))+
                               geom_histogram(position="identity",
                                              alpha=.5, bins=30)+
                               scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                               geom_vline(aes(xintercept=.y$SRMR[1],
                                              linetype="final$SRMR[1]",color="final$SRMR[1]"),
                                          size=.6)+
                               geom_vline(aes(xintercept=.08,
                                              linetype=".08",color=".08"),
                                          size=.75)+
                               scale_color_manual(name="Cutoff Values",
                                                  labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                  values=c("final$SRMR[1]"="black",
                                                           ".08"="black"))+
                               scale_linetype_manual(name="Cutoff Values",
                                                     labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                     values=c("final$SRMR[1]"="longdash",
                                                              ".08"="dotted"))+
                               theme(axis.title.y = element_blank(),
                                     axis.text.y = element_blank(),
                                     axis.ticks.y = element_blank(),
                                     panel.background = element_blank(),
                                     axis.line = element_line(color="black"),
                                     legend.position = "none",
                                     legend.title = element_blank(),
                                     legend.box = "vertical"))

    #Plot RMSEA.  Need map2 and data=.x (can't remember why).
    RMSEA_plot <- purrr::map2(plot,final,~ggplot(data=.x,aes(x=RMSEA,fill=Model))+
                                geom_histogram(position="identity",
                                               alpha=.5, bins=30)+
                                scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                                geom_vline(aes(xintercept=.y$RMSEA[1],
                                               linetype="final$RMSEA[1]",color="final$RMSEA[1]"),
                                           size=.6)+
                                geom_vline(aes(xintercept=.06,
                                               linetype=".06",color=".06"),
                                           size=.75)+
                                scale_color_manual(name="Cutoff Values",
                                                   labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                   values=c("final$RMSEA[1]"="black",
                                                            ".06"="black"))+
                                scale_linetype_manual(name="Cutoff Values",
                                                      labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                      values=c("final$RMSEA_M[1]"="longdash",
                                                               ".06"="dotted"))+
                                theme(axis.title.y = element_blank(),
                                      axis.text.y = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      panel.background = element_blank(),
                                      axis.line = element_line(color="black"),
                                      legend.position = "none",
                                      legend.title = element_blank(),
                                      legend.box = "vertical"))

    #Plot CFI. Need map2 and data=.x (can't remember why).
    CFI_plot <- purrr::map2(plot,final,~ggplot(data=.x,aes(x=CFI,fill=Model))+
                              geom_histogram(position="identity",
                                             alpha=.5, bins=30)+
                              scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                              geom_vline(aes(xintercept=.y$CFI[1],
                                             linetype="final$CFI[1]",color="final$CFI[1]"),
                                         size=.6)+
                              geom_vline(aes(xintercept=.95,
                                             linetype=".95",color=".95"),
                                         size=.75)+
                              scale_color_manual(name="Cutoff Values",
                                                 labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                 values=c("final$CFI[1]"="black",
                                                          ".95"="black"))+
                              scale_linetype_manual(name="Cutoff Values",
                                                    labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                    values=c("final$CFI[1]"="longdash",
                                                             ".95"="dotted"))+
                              theme(axis.title.y = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    panel.background = element_blank(),
                                    axis.line = element_line(color="black"),
                                    legend.position = "none",
                                    legend.title = element_blank(),
                                    legend.box = "vertical"))

    #Create a list with the plots combined for each severity level
    plots_combo <- base::lapply(base::seq(base::length(plot)),function(x) c(SRMR_plot[x],RMSEA_plot[x],CFI_plot[x]))

    #Add a collective legend and title with the level indicator
    plots <- base::lapply(base::seq(base::length(plots_combo)), function(x) patchwork::wrap_plots(plots_combo[[x]])+
                            plot_layout(guides = "collect")+
                            plot_annotation(title=paste("Level", x))
                          & theme(legend.position = 'bottom'))

    #Put into list
    res$plots <- plots

    dataMiss<-data_nnor(model9, data, n, reps)

    jj<-as.data.frame(dataMiss$data_t)[,names(dataMiss$data1)]

    aa<-gather(jj)
    bb<-gather(dataMiss$data1)

    aa$model<-c("Bootstrapped Data")
    bb$model<-c("Original Data")
    cc<-rbind(aa,bb)

    dist<- ggplot()+  geom_histogram(data=cc,aes(x=value, y=..density.., fill=model), bins=(max(bb$value, na.rm=T)-min(bb$value, na.rm=T)),alpha=.3,position="identity") + facet_wrap(~key, scales = 'free_x')+
      scale_colour_manual(values=c("#E9798C","#66C2F5")) +
      theme(axis.title.y = element_blank(), axis.title.x=element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(color="black"),legend.title=element_blank())

    res$dist_plot <- dist


  }

  #Create object (necessary for subsequent print statement)
  class(res) <- 'nnorOne'

  return(res)

}

#' @method print nnorOne
#' @param x nnorOne object
#' @param ... other print parameters
#' @rdname nnorOne
#' @export

#Print suppression/organization statement for list
#Needs same name as class, not function name
#Need to add ... param or will get error message in CMD check
print.nnorOne <- function(x,...){

  #Automatically return fit cutoffs
  base::cat("Your DFI cutoffs: \n")
  base::print(x$cutoffs)

  #Only print fit indices from lavaan object if someone submits a lavaan object
  if(!is.null(x$fit)){
    base::cat("\n")

    base::cat("Empirical fit indices: \n")
    base::print(x$fit)
  }

  base::cat("\n Notes:
  -'Sensitivity' is % of hypothetically misspecified models correctly identified by cutoff in DFI simulation
  -Cutoffs with 95% sensitivity are reported when possible
  -If sensitivity is <50%, cutoffs will be supressed \n")

  if(!is.null(x$plots)){

    base::cat("\n The distributions for each level are in the Plots tab \n")
    base::print(x$plots)
  }

  #Hides this function
  base::invisible()
}
