#' @title Dynamic fit index (DFI) cutoffs for categorical one-factor CFA models
#'
#' @description This function generates DFI cutoffs for one-factor CFA models that treat items as categorical.
#' The default argument is a singular argument: a \code{\link{lavaan}} object from the \code{\link{cfa}} function.
#' The function can also accommodate manual entry of the model statement and sample size (including threshold estimates).
#'
#' The app-based version of this function can be found at \href{https://dynamicfit.app/}{dynamicfit.app}.
#'
#' @param model This can either be a \code{\link{lavaan}} object from the \code{\link{cfa}} function,
#' OR a model statement written in \code{\link{lavaan}} \code{\link{model.syntax}} with standardized loadings and thresholds.
#' @param n If you entered a \code{\link{lavaan}} object for model, leave this blank.
#' Otherwise, enter your sample size (numeric).
#' @param plot Displays distributions of fit indices for each level of misspecification.
#' @param manual If you entered a \code{\link{lavaan}} object, keep this set to FALSE.
#' If you manually entered standardized loadings and sample size, set this to TRUE.
#' @param reps The number of replications used in your simulation. This is set to 500 by default in both the
#' R package and the corresponding Shiny App.
#' @param estimator Which estimator to use within the simulations (enter in quotes). The default is WLSMV. Only limited-information estimators that produce fit indices are permitted (i.e., maximum likelihood is not available)
#'
#' @import dplyr lavaan simstandard ggplot2 stringr rlang
#' @importFrom purrr map map_dfr map2
#' @importFrom tidyr unite extract
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#'
#' @author Daniel McNeish & Melissa G Wolf
#'
#' Maintainer: Daniel McNeish <dmcneish@asu.edu>
#'
#' @rdname catOne
#'
#' @return Dynamic fit index (DFI) cutoffs for SRMR, RMSEA, and CFI.
#' @export
#'
#' @examples
#' #Example using a lavaan object as input (manual=FALSE)
#'
#' #one-factor model
#' m1<-"F1=~X5+ X6 + X7 + X8 + X9"
#'
#'  #fit the model in lavaan, treating items are categorical
#'  fit<-lavaan::cfa(m1, data=Example, ordered=TRUE)
#'
#' \donttest{catOne(fit)}
#'
#' #Manual entry example (manual=TRUE)
#'
#' #one-factor model with correlated factors
#' m1<-"F1=~X5+ X6 + X7 + X8 + X9"
#'
#'  #fit the model, treating items are categorical
#'  #lavaan is used here to shown where estimates come from
#'  #but manual entry supports standardized estimates from models fit in any software
#'
#'  fit<-lavaan::cfa(m1, data=Example, ordered=TRUE)
#'  lavaan::standardizedsolution(fit)
#'
#' #thresholds go in model statement as
#'   #(a)categorical item name
#'   #(b) vertical pipe
#'   #(c) estimate
#'   #(d)times t+threshold number
#'
#' manual_model <-"F1=~.550*X5 + .614*X6 + .726*X7 + .723*X8 + .236*X9
#'
#'X5 |-0.274*t1
#'X5 | 0.305*t2
#'X5 | 0.765*t3
#'X5 | 1.259*t4
#'
#'X6 |-0.279*t1
#'X6 | 0.353*t2
#'X6 | 0.779*t3
#'X6 | 1.175*t4
#'
#'X7 |-0.269*t1
#'X7 | 0.385*t2
#'X7 | 0.871*t3
#'X7 | 1.329*t4
#'
#'X8 |-0.274*t1
#'X8 | 0.358*t2
#'X8 | 0.779*t3
#'X8 | 1.237*t4
#'
#'X9 |-0.269*t1
#'X9 | 0.342*t2
#'X9 | 0.745*t3
#'X9 | 1.248*t4"
#'
#' \donttest{catOne(model=manual_model,n=500,manual=TRUE)}
#'

catOne <- function(model,n=NULL,plot=FALSE,manual=FALSE,reps=250, estimator="WLSMV"){
  #If manual, expect manual (a la Shiny app)
  if(manual){

    tryCatch(modelWithNum(model),
             error=function(err5){
               if (grepl("no method for coercing this S4 class to a vector", err5)){
                 stop("dynamic Error: Did you accidentally include 'manual=TRUE' with a non-manually entered lavaan object?")
               }
             })

    n <- n
    model9 <- modelWithNum(model) #DELETE THRESHOLDS FROM MODEL STATEMENT
    threshold<-cleanthreshold(model) #COLLECT THRESHOLDS IN SEPARATE OBJECT TO PASS TO LATER FUNCTIONS

    tryCatch(defre(model9,n),
             error=function(err4){
               if (grepl("non-numeric matrix extent", err4)){
                 stop("dynamic Error: Did you forget to include a sample size with your manually entered model?")
               }
             })

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
    threshold <-Thresh(model) # ISOLATE THRESHOLDS FROM LAVAAN OBJECT TO PASS THEM TO LATER FUNCTIONS
  }

  if (unstandardized(model9)>0){
    stop("dynamic Error: One of your loadings or correlations has an absolute value of 1 or above (an inadmissible value). Please use standardized loadings. If all of your loadings are under 1, try looking for a missing decimal somewhere in your model statement.")
  }

  if (number_factor(model9)>1){
    stop("dynamic Error: You entered a multi-factor model.  Use catHB instead.")
  }

  if (defre(model9,n)==0){
    stop("dynamic Error: Misspecifications cannot be added to a just identified model.")
  }

  if ( nrow(one_num(model9)) < (number_factor(model9)-1)){
    stop("dynamic Error: There are not enough free items to produce all misspecification levels.")
  }

  if (nrow(threshold) <1){
    stop("dynamic Error: No thresholds were found. Make sure you treated data as categorical with an 'ordered=' option in lavaan or that you included threshold estimates in your model statement.")
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
  results <- one_df_cat(model9,n,reps, threshold,estimator)

  #Save the data and make it exportable
  res$data <- fit_data(results)

  #For each list element (misspecification) compute the cutoffs
  misspec_sum <- purrr::map(results,~dplyr::reframe(.,SRMR_M=stats::quantile(SRMR_M, c(seq(0.05,1,0.01))),
                                                      RMSEA_M=stats::quantile(RMSEA_M, c(seq(0.05,1,0.01))),
                                                      CFI_M=stats::quantile(CFI_M, c(seq(0.95,0,-0.01)))))

  #For the true model, compute the cutoffs (these will all be the same - just need in list form)
  true_sum <- purrr::map(results,~dplyr::reframe(.,SRMR_T=stats::quantile(SRMR_T, c(.95)),
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
    S[[i]]<-subset(fit[[i]], subset=(!duplicated(fit[[i]][('S')])|fit[[i]][('Power')]==0), select=c("SRMR_M","Power","S")) %>% filter(S==1|Power==0)
    R[[i]]<-subset(fit[[i]], subset=(!duplicated(fit[[i]][('R')])|fit[[i]][('Power')]==0), select=c("RMSEA_M","Power","R")) %>% filter(R==1|Power==0)
    C[[i]]<-subset(fit[[i]], subset=(!duplicated(fit[[i]][('C')])|fit[[i]][('Power')]==0), select=c("CFI_M","Power","C"))  %>% filter(C==1|Power==0)
    colnames(S[[i]])<-c("SRMR","PowerS")
    colnames(R[[i]])<-c("RMSEA","PowerR")
    colnames(C[[i]])<-c("CFI","PowerC")

    final[[i]]<-cbind(S[[i]][1,],R[[i]][1,],C[[i]][1,])
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
    #                                                   RMSEA_M=stats::quantile(RMSEA_M, c(.05,.1)),
    #                                                  CFI_M=stats::quantile(CFI_M, c(.95,.9))))

    #For the true model, compute the cutoffs (these will all be the same - just need in list form)
    #true_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_T=stats::quantile(SRMR_T, c(.95,.9)),
    #                                                RMSEA_T=stats::quantile(RMSEA_T, c(.95,.9)),
    #                                               CFI_T=stats::quantile(CFI_T, c(.05,.1))))

    #Select just those variables and rename columns to be the same
    Misspec_dat <- purrr::map(results,~dplyr::select(.,SRMR_M:Type_M) %>%
                                `colnames<-`(c("SRMR","RMSEA","CFI","Model")))

    #Select just those variables and rename columns to be the same
    True_dat <- purrr::map(results,~dplyr::select(.,SRMR_T:Type_T) %>%
                             `colnames<-`(c("SRMR","RMSEA","CFI", "Model")))

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
                                                  labels=c("Hu & Bentler Cutoff","Dynamic Cutoff"),
                                                  values=c("final$SRMR[1]"="black",
                                                           ".08"="black"))+
                               scale_linetype_manual(name="Cutoff Values",
                                                     labels=c("Hu & Bentler Cutoff","Dynamic Cutoff"),
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
                                                   labels=c("Hu & Bentler Cutoff","Dynamic Cutoff"),
                                                   values=c("final$RMSEA[1]"="black",
                                                            ".06"="black"))+
                                scale_linetype_manual(name="Cutoff Values",
                                                      labels=c("Hu & Bentler Cutoff","Dynamic Cutoff"),
                                                      values=c("final$RMSEA[1]"="longdash",
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
                                                 labels=c("Hu & Bentler Cutoff","Dynamic Cutoff"),
                                                 values=c("final$CFI[1]"="black",
                                                          ".95"="black"))+
                              scale_linetype_manual(name="Cutoff Values",
                                                    labels=c("Hu & Bentler Cutoff","Dynamic Cutoff"),
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

  }

  #Create object (necessary for subsequent print statement)
  class(res) <- 'catOne'

  return(res)

}


#' @method print catOne
#' @param x catOne object
#' @param ... other print parameters
#' @rdname catOne
#' @export

#Print suppression/organization statement for list
#Needs same name as class, not function name
#Need to add ... param or will get error message in CMD check
print.catOne <- function(x,...){

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
