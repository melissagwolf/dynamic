#' @title Dynamic fit index (DFI) cutoffs for one-factor CFA models
#'
#' @description This function generates DFI cutoffs for one-factor CFA models using ML estimation.
#' The default argument is a singular argument: a \code{\link{lavaan}} object from the \code{\link{cfa}} function.
#' The function can also accommodate manual entry of the model statement and sample size.
#'
#' The app-based version of this function can be found at \href{https://dynamicfit.app/}{dynamicfit.app}.
#'
#' @param model This can either be a \code{\link{lavaan}} object from the \code{\link{cfa}} function,
#' OR a model statement written in \code{\link{lavaan}} \code{\link{model.syntax}} with standardized loadings.
#' @param n If you entered a \code{\link{lavaan}} object for model, leave this blank.
#' Otherwise, enter your sample size (numeric).
#' @param plot Displays distributions of fit indices for each level of misspecification.
#' @param manual If you entered a \code{\link{lavaan}} object, keep this set to FALSE.
#' If you manually entered standardized loadings and sample size, set this to TRUE.
#' @param reps (**Do not modify this**): The number of replications used in your simulation. This is set to 500 by default in both the
#' R package and the corresponding Shiny App.
#' @param estimator Which estimator to use within the simulations (enter in quotes). The default is ML.
#'
#' @import dplyr lavaan simstandard ggplot2 stringr
#' @importFrom purrr map map_dfr map2
#' @importFrom tidyr unite extract
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#'
#' @author Melissa G Wolf & Daniel McNeish
#'
#' Maintainer: Melissa G Wolf <melissagordon@ucsb.edu>
#'
#' @rdname cfaOne
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
#' #Manual entry example for a sample size of 300 (manual=TRUE)
#' manmod <- "F1 =~ .602*Y1 + .805*Y2 + .857*Y3 + .631*Y4 + .345*Y5 + .646*Y6"
#' \donttest{cfaOne(model=manmod,n=300,manual=TRUE)}
#'

cfaOne <- function(model,n=NULL,plot=FALSE,manual=FALSE,estimator="ML",reps=500){

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
    fitted <- round(lavaan::fitmeasures(model,c("chisq","df","pvalue","srmr","rmsea","cfi")),3)
    fitted_m <- as.matrix(fitted)
    fitted_t <- t(fitted_m)
    fitted_t <- as.data.frame(fitted_t)
    colnames(fitted_t) <- c("Chi-Square"," df","p-value","  SRMR","  RMSEA","   CFI")
    rownames(fitted_t) <- c("")
    res$fit <- fitted_t
  }

  #Run simulation
  results <- one_df(model9,n,reps)

  #Save the data and make it exportable
  res$data <- fit_data(results)

  #For each list element (misspecification) compute the cutoffs
  misspec_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_M=stats::quantile(SRMR_M, c(.05,.1)),
                                                      RMSEA_M=stats::quantile(RMSEA_M, c(.05,.1)),
                                                      CFI_M=stats::quantile(CFI_M, c(.95,.9))))

  #For the true model, compute the cutoffs (these will all be the same - just need in list form)
  true_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_T=stats::quantile(SRMR_T, c(.95,.9)),
                                                   RMSEA_T=stats::quantile(RMSEA_T, c(.95,.9)),
                                                   CFI_T=stats::quantile(CFI_T, c(.05,.1))))

  #Bind each of the misspecified cutoffs to the true cutoffs, listwise
  Table <- purrr::map(misspec_sum,~base::cbind(.,true_sum[[1]]) %>%
                        dplyr::mutate(SRMR_R=base::round(SRMR_M,3),
                                      RMSEA_R=base::round(RMSEA_M,3),
                                      CFI_R=base::round(CFI_M,3),
                                      SRMR=base::ifelse(SRMR_T<SRMR_M,SRMR_R,"NONE"),
                                      RMSEA=base::ifelse(RMSEA_T<RMSEA_M,RMSEA_R,"NONE"),
                                      CFI=base::ifelse(CFI_T>CFI_M,CFI_R,"NONE")) %>%
                        dplyr::select(SRMR,RMSEA,CFI))

  #This is to clean up the table for presentation
  #list is a function within mutate to apply function lead across each element
  Row2 <- purrr::map_dfr(Table,~dplyr::mutate(.,SRMR_1=SRMR,
                                              RMSEA_1=RMSEA,
                                              CFI_1=CFI) %>%
                           dplyr::mutate_at(c("SRMR_1","RMSEA_1","CFI_1"),base::list(dplyr::lead)) %>%
                           dplyr::slice(1) %>%
                           dplyr::mutate(SRMR=base::ifelse(base::is.character(SRMR),SRMR_1,"--"),
                                         RMSEA=base::ifelse(base::is.character(RMSEA),RMSEA_1,"--"),
                                         CFI=base::ifelse(base::is.character(CFI),CFI_1,"--"),
                                         SRMR=stringr::str_replace_all(base::as.character(SRMR),"0\\.","."),
                                         RMSEA=stringr::str_replace_all(base::as.character(RMSEA),"0\\.","."),
                                         CFI=stringr::str_replace_all(base::as.character(CFI),"0\\.",".")) %>%
                           dplyr::select(SRMR,RMSEA,CFI))

  #Still cleaning
  #Unlist Table
  Table_C <- purrr::map_dfr(Table,~dplyr::mutate(.,SRMR=stringr::str_replace_all(base::as.character(SRMR),"0\\.","."),
                                                 RMSEA=stringr::str_replace_all(base::as.character(RMSEA),"0\\.","."),
                                                 CFI=stringr::str_replace_all(base::as.character(CFI),"0\\.",".")))

  #Cleaning
  Table_C[base::seq(2,nrow(Table_C),by=2),] <- Row2

  #Create row names for level
  Table_C$levelnum <- base::paste("Level", base::rep(1:(base::nrow(Table_C)/2),each=2))

  #Create row names for proportions
  Table_C$cut <- base::rep(c("95/5","90/10"))

  #Add rownames to final table
  Final_Table <- Table_C %>%
    tidyr::unite(Cut,levelnum,cut,sep=": ") %>%
    tibble::column_to_rownames(var='Cut')

  #Put into list
  res$cutoffs <- Final_Table

  #If user selects plot = T
  if(plot){

    #For each list element (misspecification) compute the cutoffs
    misspec_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_M=stats::quantile(SRMR_M, c(.05,.1)),
                                                        RMSEA_M=stats::quantile(RMSEA_M, c(.05,.1)),
                                                        CFI_M=stats::quantile(CFI_M, c(.95,.9))))

    #For the true model, compute the cutoffs (these will all be the same - just need in list form)
    true_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_T=stats::quantile(SRMR_T, c(.95,.9)),
                                                     RMSEA_T=stats::quantile(RMSEA_T, c(.95,.9)),
                                                     CFI_T=stats::quantile(CFI_T, c(.05,.1))))

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
    SRMR_plot <- purrr::map2(plot,misspec_sum,~ggplot(data=.x,aes(x=SRMR,fill=Model))+
                               geom_histogram(position="identity",
                                              alpha=.5, bins=30)+
                               scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                               geom_vline(aes(xintercept=.y$SRMR_M[1],
                                              linetype="misspec_sum$SRMR_M[1]",color="misspec_sum$SRMR_M[1]"),
                                          size=.6)+
                               geom_vline(aes(xintercept=.08,
                                              linetype=".08",color=".08"),
                                          size=.75)+
                               scale_color_manual(name="Cutoff Values",
                                                  labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                  values=c("misspec_sum$SRMR_M[1]"="black",
                                                           ".08"="black"))+
                               scale_linetype_manual(name="Cutoff Values",
                                                     labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                     values=c("misspec_sum$SRMR_M[1]"="longdash",
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
    RMSEA_plot <- purrr::map2(plot,misspec_sum,~ggplot(data=.x,aes(x=RMSEA,fill=Model))+
                                geom_histogram(position="identity",
                                               alpha=.5, bins=30)+
                                scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                                geom_vline(aes(xintercept=.y$RMSEA_M[1],
                                               linetype="misspec_sum$RMSEA_M[1]",color="misspec_sum$RMSEA_M[1]"),
                                           size=.6)+
                                geom_vline(aes(xintercept=.06,
                                               linetype=".06",color=".06"),
                                           size=.75)+
                                scale_color_manual(name="Cutoff Values",
                                                   labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                   values=c("misspec_sum$RMSEA_M[1]"="black",
                                                            ".06"="black"))+
                                scale_linetype_manual(name="Cutoff Values",
                                                      labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                      values=c("misspec_sum$RMSEA_M[1]"="longdash",
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
    CFI_plot <- purrr::map2(plot,misspec_sum,~ggplot(data=.x,aes(x=CFI,fill=Model))+
                              geom_histogram(position="identity",
                                             alpha=.5, bins=30)+
                              scale_fill_manual(values=c("#E9798C","#66C2F5"))+
                              geom_vline(aes(xintercept=.y$CFI_M[1],
                                             linetype="misspec_sum$CFI_M[1]",color="misspec_sum$CFI_M[1]"),
                                         size=.6)+
                              geom_vline(aes(xintercept=.95,
                                             linetype=".95",color=".95"),
                                         size=.75)+
                              scale_color_manual(name="Cutoff Values",
                                                 labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                 values=c("misspec_sum$CFI_M[1]"="black",
                                                          ".95"="black"))+
                              scale_linetype_manual(name="Cutoff Values",
                                                    labels=c("Hu & Benter Cutoff","Dynamic Cutoff"),
                                                    values=c("misspec_sum$CFI_M[1]"="longdash",
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
  class(res) <- 'cfaOne'

  return(res)

}

#' @method print cfaOne
#' @param x cfaOne object
#' @param ... other print parameters
#' @rdname cfaOne
#' @export

#Print suppression/organization statement for list
#Needs same name as class, not function name
#Need to add ... param or will get error message in CMD check
print.cfaOne <- function(x,...){

  #Automatically return fit cutoffs
  base::cat("Your DFI cutoffs: \n")
  base::print(x$cutoffs)

  #Only print fit indices from lavaan object if someone submits a lavaan object
  if(!is.null(x$fit)){
    base::cat("\n")

    base::cat("Empirical fit indices: \n")
    base::print(x$fit)
  }

  if(!is.null(x$plots)){

    base::cat("\n The distributions for each level are in the Plots tab \n")
    base::print(x$plots)
  }

  #Hides this function
  base::invisible()
}
