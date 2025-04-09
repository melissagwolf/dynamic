#' @title Dynamic fit index (DFI) cutoffs multi-factor CFA models with Likert-type items
#'
#' @description This function generates DFI cutoffs for multi-factor CFA models that treats Likert-type items as continuous.
#'  The default argument is a singular argument: a \code{\link{lavaan}} object from
#'  the \code{\link{cfa}} function. The function can also accommodate manual entry of the model statement and
#'  sample size. A primary difference in likert DFI functions is that a dataset must also be provided in the 'data' argument
#'  in order to simulate data with the same number of response options and response frequencies as the original data.
#'  When Likert-type items are treated as continuous, this information cannot be obtained solely from model output.
#'
#' The app-based version of this function can be found at \href{https://dynamicfit.app/}{dynamicfit.app}.
#'
#' @param model This can either be a \code{\link{lavaan}} object from the \code{\link{cfa}} function,
#' OR a model statement written in \code{\link{lavaan}} \code{\link{model.syntax}} with standardized loadings.
#' @param data An empirical dataset used to determine the number of Likert responses and the response frequencies
#' @param n If you entered a \code{\link{lavaan}} object for model, leave this blank.
#' Otherwise, enter your sample size (numeric).
#' @param plot Displays distributions of fit indices for each level of misspecification.
#' This also includes plots to visualize how close the simulated data are to the original data.
#' @param manual If you entered a \code{\link{lavaan}} object, keep this set to FALSE (the default).
#' If you manually entered standardized loadings and sample size, set this to TRUE.
#' @param reps The number of replications used in your simulation. This is set to 500 by default.
#' @param estimator Which estimator to use within the simulations (enter in quotes). The default is ML
#'
#' @import dplyr lavaan simstandard ggplot2 stringr
#' @importFrom purrr map map_dfr map2
#' @importFrom tidyr unite extract gather
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#'
#' @author Daniel McNeish & Melissa G Wolf
#'
#' Maintainer: Daniel McNeish <dmcneish@asu.edu>
#'
#' @rdname likertHB
#'
#' @return Dynamic fit index (DFI) cutoffs for SRMR, RMSEA, and CFI.
#' @export

#' @examples
#' #Example using a lavaan object as input (manual=FALSE)
#'
#' #two-factor model with correlated factors
#' m1<-"
#'  F1=~X2 + X3 + X4
#'  F2=~X6 + X7 + X8
#'  F1~~F2"
#'
#'  #fit the model in lavaan, treating items are continuous
#'  fit<-lavaan::cfa(m1, data=Example)
#'
#' \donttest{likertHB(fit, data=Example)}
#'
#' #Manual entry example (manual=TRUE)
#'
#' #two-factor model with correlated factors
#' m1<-"
#'  F1=~X2 + X3 + X4
#'  F2=~X6 + X7 + X8
#'  F1~~F2"
#'
#'  #fit the model, treating items are continuous
#'  #lavaan is used here to shown where estimates come from
#'  #but manual entry supports standardized estimates from models fit in any software
#'
#'  fit<-lavaan::cfa(m1, data=Example)
#'  lavaan::standardizedsolution(fit)
#'
#'
#' manual_model <-"F1=~.554*X2 + .654*X3 + .733*X4
#'  F2=~.537*X6 + .666*X7 + .723*X8
#'  F1~~.339*F2"
#'
#' \donttest{likertHB(model=manual_model,data=Example,n=500,manual=TRUE)}
#'

likertHB <- function(model,data,n=NULL,plot=FALSE,manual=FALSE,estimator="ML",reps=250){

  #If manual, expect manual (a la Shiny app)
  if(manual){

    tryCatch(cleanmodel(model),
             error=function(err5){
               if (grepl("no method for coercing this S4 class to a vector", err5)){
                 stop("dynamic Error: Did you accidentally include 'manual=TRUE' with a non-manually entered lavaan object?")
               }
             })

    n <- n
    model9 <- model

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
  }

  if (unstandardized(model9)>0){
    stop("dynamic Error: One of your loadings or correlations has an absolute value of 1 or above (an inadmissible value). Please use standardized loadings. If all of your loadings are under 1, try looking for a missing decimal somewhere in your model statement.")
  }

  if (number_factor(model9)<2){
    stop("dynamic Error: You entered a one-factor model.  Use likertOne function instead.")
  }


  if (defre(model9,n)==0){
    stop("dynamic Error: Misspecifications cannot be added to a just identified model.")
  }

  if ( nrow(multi_num_HB(model9)) < (number_factor(model9)-1)){
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
  tryCatch(multi_df_likert(model9,data,n,estimator, reps),
           error=function(err3){
             if (grepl("is missing, with no default", err3)){
               stop("dynamic Error: Did you forget to include a dataset? likert functions in dynamic require the original data to accurately simulate data that have the same
                      response proportions as the original data. If you do not have the data, you can use the cfaHB function,
                      although this will assume truly continuous responses rather than Likert responses. ")
             }
           })
  results <- multi_df_likert(model9,data,n,estimator, reps)

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

  ##append table with other cutoffs,  loop to control different number of levels
  Fit<-round(L0,3)
  for (i in 1:length(misspec_sum)){
    Fit<-round(rbind(Fit, final[[i]]),3)
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

  #create blanks to make table easier to read
  pp<-c(rep("--",(length(misspec_sum)+1)))
  pp0<-c(rep("",(length(misspec_sum)+1)))

  #matrix of cross-loadings added at each level
  mag <- multi_add_HB(model9) %>%
    tidyr::separate(V1,into=c("a","b","Magnitude","d","e"),sep=" ") %>%
    select(Magnitude) %>%
    mutate(Magnitude=as.numeric(Magnitude),
           Magnitude=round(Magnitude,digits=3))

  #Create column of cross-loadings added at each level
  mname<-c("NONE","","")
  for (i in 1:length(misspec_sum)){
    mi<-c(mag[i,],"","")
    mname<-cbind(mname,mi)
  }

  #format column for each index and the misspecification magnitude
  SS<-noquote(matrix(rbind(fit1[,1],fit1[,2],pp0),ncol=1))
  RR<-noquote(matrix(rbind(fit1[,3],fit1[,4],pp0),ncol=1))
  CC<-noquote(matrix(rbind(fit1[,5],fit1[,6],pp0),ncol=1))
  MM<-noquote(matrix(mname, ncol=1))

  #bind columns together into one table
  Table<-noquote(cbind(SS,RR,CC,MM) %>%
                   `colnames<-`(c("SRMR","RMSEA","CFI","Magnitude")))

  #update row names
  rname<-c("Level-0","Specificity", "")
  for (i in 1:length(misspec_sum)){
    ri<-c(paste("Level-",i,sep=""),"Sensitivity","")
    rname<-cbind(rname,ri)
  }
  rownames(Table)<-rname

  #delete last blank row
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

    dataMiss<-data_likert(model9, data, n)



    aa<-tidyr::gather(dataMiss$sim)
    bb<-tidyr::gather(dataMiss$orig)

    aa$model<-c("Simulated Data")
    bb$model<-c("Original Data")
    cc<-rbind(aa,bb)

    dist<- ggplot()+  geom_bar(data=cc,aes(x=value, fill=model),alpha=.3,position="identity") + facet_wrap(~key, scales = 'free_x')+
      scale_colour_manual(values=c("#E9798C","#66C2F5")) + scale_x_continuous(breaks=min(bb$value, na.rm=T):(max(bb$value, na.rm=T)))+
      theme(axis.title.y = element_blank(), axis.title.x=element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(color="black"),legend.title=element_blank())

    res$dist_plot <- dist


  }

  #Create object (necessary for subsequent print statement)
  class(res) <- 'likertHB'

  return(res)

}

#' @method print likertHB
#' @param x likertHB object
#' @param ... other print parameters
#' @rdname likertHB
#' @export

#Print suppression/organization statement for list
#Needs same name as class, not function name
#Need to add ... param or will get error message in CMD check
print.likertHB <- function(x,...){

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
