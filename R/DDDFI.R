#' @title Direct Discrepancy Dynamic fit index (3DFI) cutoffs for arbitrary covariance structure models
#'
#' @description This function generates DFI cutoffs for any single group covariance structure model with
#' a saturated or absent mean structure. It supports (a) any estimator supported by lavaan (e.g.
#' ML, MLR, WLSMV, ULSMV), (b) missing data, and (c) multiple response scales (normal, non-normal continuous,
#' categorical). The default argument is a singular argument: a \code{\link{lavaan}} object.The function
#' can also accommodate manual entry of the model statement and sample size. Some features require an original
#' dataset to be provided (e.g., missing data, categorical data). The app-based version of
#' this function can be found at \href{https://dynamicfit.app/}{dynamicfit.app}.
#'
#' @param model This can either be a \code{\link{lavaan}} object, OR a model statement written in \code{\link{lavaan}} \code{\link{model.syntax}}
#' with standardized estimates
#' @param scale Determines how data are simulated. Options are "normal", "nonnormal", or "categorical". "normal" assumes multivariate
#' normality across all variables. "nonnormal" recreates distributions in an empirical dataset (to be provided by the user), assuming variables are continuous.
#' "categorical" simulates discrete data with the same proportions as an empirical dataset (to be provided by the user). With "categorical", mixed formats
#' are also supported and any variable with more than 9 categories is simulated from a normal distribution.Only "normal" can be used without provided an original
#' dataset.
#' @param data The original data to which the model was applied. Not required if scale="normal". Otherwise, data is required.
#' @param manual If you entered a \code{\link{lavaan}} object, keep this set to FALSE.
#' If you manually entered standardized estimates and sample size, set this to TRUE.
#' @param reps The number of replications used in the simulations. This is set to 250 by default
#' @param n If you entered a \code{\link{lavaan}} object for model, leave this blank.
#' Otherwise, enter your sample size (numeric).
#' @param estimator Which estimator to use within the simulations (enter in quotes). The default depends on the scale option ("ML" for "normal", "MLR" for "nonnormal",
#' and "WLSMV" for categorical)
#' @param MAD Mean Absolute Discrepancies to test in the simulation. Default is c(.038, .05, .06) to recreate traditional "Close", "Fair", "Mediocre" benchmarks
#' @param plot.dfi Displays simulated distributions of fit indices used to derive cutoffs for each MAD value.
#' @param plot.dist Displays distributions of simulated data (and empirical data, if provided) to assess fidelity of simulated data to empirical data
#' @param plot.discrepancy Displays distributions of simulated MAD values
#'
#' @import dplyr lavaan simstandard ggplot2 stringr
#' @importFrom purrr map map_dfr map2
#' @importFrom tidyr unite extract
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#' @importFrom semTools bsBootMiss
#' @importFrom MASS mvrnorm

#'
#' @author Daniel McNeish & Melissa G Wolf
#'
#' Maintainer: Daniel McNeish <dmcneish@asu.edu>
#'
#' @rdname DDDFI
#'
#' @return Direct Discrepancy Dynamic fit index (DFI) cutoffs for CFI, RMSEA, and RMSEA 90% CI.
#' @export
#'


DDDFI <- function(model,data=NULL, scale="normal", manual=FALSE,reps=250, n=NULL, estimator=NULL, MAD=c(.038, .05,.06), plot.dfi=FALSE, plot.dist=FALSE, plot.discrepancy=FALSE){
if (scale=="normal" & is.null(estimator)){
  estimator="ML"
}

if (scale=="nonnormal" & is.null(estimator)){
  estimator="MLR"
}

if (scale=="categorical" & is.null(estimator)){
  estimator="WLSMV"
}

if (scale %in% c("nonnormal", "categorical") & is.null(data)){
    stop("dynamic error: A dataset is required for scale='nonnormal' or scale='categorical'

         If you do not have a dataset, scale='normal' may be used without a dataset.

         Without a dataset, cutoffs may be affected by missing data and deviations from normality.")
  }

if ( !(scale%in% c("normal", "nonnormal", "categorical"))){
  stop("dynamic Error:  scale must be one of 'normal','nonnormal', or 'categorical'")}

if(plot.dist & is.null(data)){
    stop("dynamic Error: Plots comparing generated data to the original data require a dataset to be provided.

           Please provide a original data with the 'data=' option or set the 'plot.dist=' option to FALSE
           ")
  }

#If model statement is manually entered,
if(manual){

  # error for using manual=T with a lavaan object input
  tryCatch(cleanmodel_DDFI(model),
           error=function(err5){
             if (grepl("no method for coercing this S4 class to a vector", err5)){
               stop("dynamic Error: Did you accidentally include 'manual=TRUE' with a non-manually entered lavaan object?")
             }
           })

  #error for forgetting to include a sample size with manual=T & scale=Normal
  if (scale=="normal" & is.null(n)){
    stop("dynamic Error: A sample size is required in the n= option with manual = T and scale='normal'")
  }

  if (!is.null(n)){
    n <- n
  }

  if (scale %in% c("nonnormal", "categorical") & is.null(n)){
    n <- nrow(data)
  }
  model9 <- model
}

#for model statement entered from a lavaan object,
if(!manual){

  #error for entering a model statement manually but forgetting to include manual=T
  if (is.character(model)){
    stop("dynamic error: Did you forget to use manual=TRUE?")
  }

  #Error for when someone enters an object that doesn't exist, or a non-lavaan object
  tryCatch(cfa_n(model),
           error=function(err2){
             if (grepl("Error in base::unlist", err2)){
               stop("dynamic Error: Did you enter a lavaan object? Confirm that it is a lavaan object using class(). If you do not have a lavaan object, enter the arguments manually and select manual=TRUE.")
             }
           })

  if (scale=="normal"){
    n <-unlist(model@Data@nobs)
  }

  if (!(scale %in% c("normal"))) {
    n<-unlist(model@Data@norig)
  }

  #convert lavaan object to manual statement upon which functions are based
  model9 <- cfa_lavmod(model)
}

####
#general error regardless of input method
####

if (unstandardized(model9)>0){
  stop("dynamic Error: One of your loadings or correlations has an absolute value of 1 or above (an inadmissible value). Please use standardized loadings. If all of your loadings are under 1, try looking for a missing decimal in your model statement or checking if you received a warning about non-positive definiteness.")
}

if (scale %in% "nonnormal"){
    warning("dynamic Warning:

            Computational times are longer for scale='nonnormal' due to missing data and procedures to generate non-normal data.

            This may take a few minutes.", immediate.=T)
  }

if (scale %in% "categorical"){
    warning("dynamic Warning:

            Computational times are longer for scale='categorical' due increased demand for categorical models.

            This may take a few minutes.", immediate.=T)
  }

#Create list to store outputs (table and plot)
res <- list()

#if input is lavaan object, extract the fit indices for the fitted model
#lots of code here because fit indices are stored in a different location/under a different name depending on the estimator/correction
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
W <- combined(model9,data,n,reps,estimator,MAD, scale)
results<-W$Table
QQ<-W$Q

results<- rapply(results, f=function(x) ifelse(is.na(x),0,x), how="replace")

#Save the data and make it exportable
res$indices <- fit_data(results)

#For each list element (misspecification) compute the cutoffs
misspec_sum <- purrr::map(results,~dplyr::reframe(.,RMSEA_CI_UPPER_M=stats::quantile(RMSEA_CI_UPPER_M, c(seq(0.05,1,0.01))),
                                                  RMSEA_M=stats::quantile(RMSEA_M, c(seq(0.05,1,0.01))),
                                                  CFI_M=stats::quantile(CFI_M, c(seq(0.95,0,-0.01)))))

#For the true model, compute the cutoffs (these will all be the same - just need in list form)
true_sum <- purrr::map(results,~dplyr::reframe(.,RMSEA_CI_UPPER_T=stats::quantile(RMSEA_CI_UPPER_T, c(.95)),
                                               RMSEA_T=stats::quantile(RMSEA_T, c(.95)),
                                               CFI_T=stats::quantile(CFI_T, c(.05))))

fit<-list()
RCI<-list()
R<-list()
C<-list()
final<-list()
for (i in 1:length(misspec_sum))
{
  fit[[i]]<-cbind(misspec_sum[[i]], true_sum[[i]])
  fit[[i]]$Power<-seq(.95, 0.0, -.01)
  fit[[i]]$R<-ifelse(fit[[i]]$RMSEA_M >= fit[[i]]$RMSEA_T, 1, 0)
  fit[[i]]$C<-ifelse(fit[[i]]$CFI_M <= fit[[i]]$CFI_T, 1, 0)

  RCI[[i]]<-subset(fit[[i]], subset=!duplicated(fit[[i]][('R')]), select=c("RMSEA_CI_UPPER_M","Power","R")) %>% filter(R==1)
  R[[i]]<-subset(fit[[i]], subset=!duplicated(fit[[i]][('R')]), select=c("RMSEA_M","Power","R")) %>% filter(R==1)
  C[[i]]<-subset(fit[[i]], subset=!duplicated(fit[[i]][('C')]), select=c("CFI_M","Power","C"))  %>% filter(C==1)

  colnames(R[[i]])<-c("RMSEA","PowerR")
  colnames(C[[i]])<-c("CFI","PowerC")
  colnames(RCI[[i]])<-c("RCI","power")

  final[[i]]<-cbind(R[[i]],C[[i]],RCI[[i]])
  final[[i]]<-final[[i]][c("RMSEA","PowerR","CFI","PowerC","RCI")]
}


L0<-data.frame(cbind(true_sum[[1]]$RMSEA_T,0.95,true_sum[[1]]$CFI_T,0.95,true_sum[[1]]$RMSEA_CI_UPPER_T))%>%
  `colnames<-`(c("RMSEA","PowerR","CFI","PowerC","RCI"))

##append table with other cutoffs,  loop to control different number of levels
Fit<-round(L0,3)
for (i in 1:length(misspec_sum)){
  Fit<-round(rbind(Fit, final[[i]]),3)
}

fit1<-unlist(Fit)%>% matrix(nrow=(length(misspec_sum)+1), ncol=5) %>%
  `colnames<-`(c("RMSEA","PowerR","CFI","PowerC","RCI"))

PR<-paste(round(100*fit1[,2], 2), "%", sep="")
PC<-paste(round(100*fit1[,4], 2), "%", sep="")

for (j in 2:(length(misspec_sum)+1)) {
  if(fit1[j,2]<.50){fit1[j,1]<-"NONE"}
  if(fit1[j,4]<.50){fit1[j,3]<-"NONE"}
  if(fit1[j,4]<.50){fit1[j,5]<-"NONE"}
}

fit1[,2]<-PR
fit1[,4]<-PC

#create blanks to make table easier to read
pp<-c(rep("--",(length(misspec_sum)+1)))
pp0<-c(rep("",(length(misspec_sum)+1)))

#create column for desired effect size
ename<-c("0.00","","")
for (i in 1:length(misspec_sum)){
  ei<-c(round(MAD[i],3),"","")
  ename<-cbind(ename,ei)
}
QQQ<-list()
for (i in 1:length(misspec_sum)){
  QQQ[[i]]<-unlist(QQ[[i]])
  QQQ[[i]]<-QQQ[[i]][QQQ[[i]] !=0]
}
actual<-unlist(lapply(QQQ,function(x) mean(abs(x))))

aname<-c("0.00","","")
for (i in 1:length(misspec_sum)){
  ai<-c(round(actual[i],3),"","")
  aname<-cbind(aname,ai)
}


#format column for each index and the misspecification magnitude
RR<-noquote(matrix(rbind(fit1[,1],fit1[,2],pp0),ncol=1))
CC<-noquote(matrix(rbind(fit1[,3],fit1[,4],pp0),ncol=1))
II<-noquote(matrix(rbind(fit1[,5],pp0,pp0),ncol=1))
EE<-noquote(matrix(ename, ncol=1))
AA<-noquote(matrix(aname,ncol=1))

#bind columns together into one table
Table<-noquote(cbind(EE,AA,CC,RR,II) %>%
                 `colnames<-`(c("MAD","Sim. MAD","CFI","RMSEA","90% CI")))


#update row names
if(all(MAD%in%c(.038,.05,.06))){
  rname<-c("Consistent","Specificity", "",
           "Close","Sensitivity", "",
           "Fair","Sensitivity", "",
           "Mediocre","Sensitivity", "")
} else{
  rname<-c("Misspecification","Specificity", "")
  for (i in 1:length(misspec_sum)){
    ri<-c("Misspecification","Sensitivity","")
    rname<-cbind(rname,ri)
  }
}
rownames(Table)<-rname

#delete last blank row
Table<-Table[1:(nrow(Table)-1),]

#Put into list
res$cutoffs <- Table

#save all simulated data
datalist<-list()
#first list element is data when fitted model is consistent with data generation model
datalist[[1]]<-dplyr::bind_rows(W$D)

#new column indicating level of misspecification
if(all(MAD%in%c(.038,.05,.06))){
  datalist[[1]]$miss_level<-"Consistent"} else{
    datalist[[1]]$miss_level<-0}

#names for default levels
mad_level<-c("Close","Fair","Medicore")
#combine list across all replications and label misspecification level
for (dl in 1:  length(misspec_sum)){
  datalist[[dl+1]]<-dplyr::bind_rows(W$Data[[2]][[dl]])
  if(all(MAD%in%c(.038,.05,.06))){
    datalist[[dl+1]]$miss_level<-mad_level[dl]} else{
      datalist[[dl+1]]$miss_level<-dl}
}

res$simulated_data<-datalist

res$estimator<-estimator
res$n<-n

#If user selects plot = T
if(plot.dfi){

  #Select just those variables and rename columns to be the same
  Misspec_dat <- purrr::map(results,~dplyr::select(.,RMSEA_M:Type_M) %>%
                              `colnames<-`(c("RMSEA","CFI","Model")))

  #Select just those variables and rename columns to be the same
  True_dat <- purrr::map(results,~dplyr::select(.,RMSEA_T:Type_T) %>%
                           `colnames<-`(c("RMSEA","CFI","Model")))

  #For each element in the list, bind the misspecified cutoffs to the true cutoffs
  #rbind doesn't work well with lists (needs do.call statement)
  plot <- base::lapply(base::seq(base::length(Misspec_dat)),function(x) dplyr::bind_rows(Misspec_dat[x],True_dat[x]))

  #Plot RMSEA.  Need map2 and data=.x (can't remember why).
  RMSEA_plot <- purrr::map2(plot,final,~ggplot(data=.x,aes(x=RMSEA,fill=Model))+
                              geom_histogram(position="identity",
                                             alpha=.5, bins=30)+
                              scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                              geom_vline(aes(xintercept=.y$RMSEA[1],
                                             linetype="final$RMSEA[1]",color="final$RMSEA[1]"),
                                         linewidth=.6)+
                              geom_vline(aes(xintercept=.06,
                                             linetype=".06",color=".06"),
                                         linewidth=.75)+
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
                                       linewidth=.6)+
                            geom_vline(aes(xintercept=.95,
                                           linetype=".95",color=".95"),
                                       linewidth=.75)+
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
  plots_combo <- base::lapply(base::seq(base::length(plot)),function(x) c(RMSEA_plot[x],CFI_plot[x]))

  #Add a collective legend and title with the level indicator
  plots <- base::lapply(base::seq(base::length(plots_combo)), function(x) patchwork::wrap_plots(plots_combo[[x]])+
                          plot_layout(guides = "collect")+
                          (if(all(MAD%in%c(.038,.05,.06))){
                            plot_annotation(title=paste(mad_level[x]))
                          } else {
                            plot_annotation(title=paste("Level", x))
                          })
                        & theme(legend.position = 'bottom'))

  #Put into list
  res$plot.dfi <- plots
}

if(plot.dist & !is.null(data)){

  gen<-dplyr::bind_rows(W$D)

  last<-ncol(gen)
  cnames<-colnames(gen[,-last])

  gen1<-gen[,-last]
  orig<-data[,cnames]

  aa<-tidyr::gather(gen1, na.rm=T)
  bb<-tidyr::gather(orig, na.rm=T)

  aa$model<-c("Generated Data")
  bb$model<-c("Original Data")
  cc<-rbind(aa,bb)
  cc$value1<-as.numeric(cc$value)

  dist<- ggplot()+  geom_histogram(data=cc,aes(x=value1, y=after_stat(density), fill=model),alpha=.3,bins=7,position="identity") + facet_wrap(~key, scales = 'free_x')+
    scale_colour_manual(values=c("#E9798C","#66C2F5")) +
    theme(axis.title.y = element_blank(), axis.title.x=element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color="black"),legend.title=element_blank())

  res$plot.dist <- dist

  ###############################
  #xxcor[lower.tri(xxcor)]

  #dummy<-sem(mod,meanstructure=T,data=gen[,colnames(a)],missing="ML",do.fit=F)
  #em<-lavInspect(dummymod,what="sampstat")
  #x<-chol(em$cov[colnames(a),colnames(a)])
  #xx<-t(x)%*%x
  #xxcor<-cov2cor(xx)

  #hist(abs(xxcor[lower.tri(xxcor)]-a[lower.tri(a)]), breaks=seq(0,.15,.01))
  #gencor<-cor(gen[,1:13],use="complete.obs")
  #hist(abs(gencor[lower.tri(gencor)]-a[lower.tri(a)]), breaks=seq(0,.15,.02))
}

if(plot.discrepancy){

  displot<-list()
  for (i in 1:length(misspec_sum)){
    displot[[i]]<-unlist(QQ[[i]])
    displot[[i]]<-as.data.frame(displot[[i]][displot[[i]] !=0])

    displot[[i]][,2]<-i
    if(all(MAD%in%c(.038,.05,.06))){
      displot[[i]][,2]<-mad_level[i]} else{
        displot[[i]][,2]<-paste("Level",i)
      }
    colnames(displot[[i]])<-c("disc","level")
  }

  displot1<-dplyr::bind_rows(displot)

  discrep<-ggplot()+ geom_histogram(data=displot1,aes(x=abs(disc),y=after_stat(density)),bins=15,alpha=.3,position="identity")+
    xlab("Discrepancy")+ggtitle("Histogram(s) of Simulated Absolute Discrepancies")+ geom_vline(xintercept=mean(abs(displot1$disc)),linetype="dotted")+ facet_wrap(~level)+
    theme(axis.title.y = element_blank(), axis.title.x=element_text(face="bold", size=12),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color="black"),legend.title=element_blank(),panel.spacing = unit(2, "lines"))

  res$plot.discrepancy<-discrep
}

#Create object (necessary for subsequent print statement)
class(res) <- 'DDDFI'

return(res)

}

#' @method print DDDFI
#' @param x DDDFI object
#' @param ... other print parameters
#' @rdname DDDFI
#' @export

#Print suppression/organization statement for list
#Needs same name as class, not function name
#Need to add ... param or will get error message in CMD check
print.DDDFI <- function(x,...){

  #Automatically return fit cutoffs
  base::cat("Your DFI cutoffs: \n")
  base::print(x$cutoffs)
  base::cat("\n")
  base::cat("Estimator:", x$estimator, "\n")
  base::cat("Sample Size:", x$n,"\n")

  #Only print fit indices from lavaan object if someone submits a lavaan object
  if(!is.null(x$fit)){
    base::cat("\n")

    base::cat("Empirical fit indices: \n")
    base::print(x$fit)

  }

  base::cat("\n Notes:
  -'Sensitivity' is % of incorrect models identified by cutoff while rejecting <5% of correct models
  - If sensitivity is <50%, cutoffs will be supressed

  -'90% CI' column is RMSEA cutoff, including sampling variability
  - Only compare upper limit of 90% RMSEA confidence interval to '90% CI' column
  - Do not compare RMSEA point estimate to '90% CI' column

  -'MAD' is the desired mean absolute discrepancy
  -'Sim. MAD' is the MAD that was achieved in the simulations
  - Values may diff when avoiding non-positive definite matrices
   \n")

  if(!is.null(x$plots)){

    base::cat("\n The distributions for each level are in the Plots tab \n")
    base::print(x$plots)
  }

  #Hides this function
  base::invisible()
}
