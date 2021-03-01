#' @title Dynamic fit index (DFI) cutoffs adapted from Hu & Bentler (1999) for multi-factor CFA models
#'
#' @description This function generates DFI cutoffs adapted from Hu & Bentler (1999) for multi-factor CFA models.
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
#' @param string If you entered a \code{\link{lavaan}} object, keep this set to FALSE.
#' If you manually entered standardized loadings and sample size, set this to TRUE.
#'
#' @import dplyr lavaan simstandard ggplot2 stringr
#' @importFrom purrr map map_dfr map2
#' @importFrom tidyr unite separate extract
#'
#' @author Melissa G Wolf & Daniel McNeish
#'
#' Maintainer: Melissa G Wolf <melissagordon@ucsb.edu>
#'
#' @rdname cfaHB
#'
#' @return Dynamic fit index (DFI) cutoffs for SRMR, RMSEA, and CFI.
#' @export
#'
#' @examples
#' #Lavaan object example (string=FALSE)
#' dat <- lavaan::HolzingerSwineford1939
#' lavmod <- "F1 =~ x1 + x2 + x3
#'            F2 =~ x4 + x5 + x6
#'            F3 =~ x7 + x8 + x9"
#' fit <- lavaan::cfa(lavmod,dat)
#' cfaHB(fit)
#'
#' #Manual entry example (string=TRUE)
#' manmod <- "F1 =~ .602*Y1 + .805*Y2 + .516*Y3 + .857*Y4
#'            F2 =~ .413*Y5 + -.631*Y6
#'            F1 ~~ .443*F2
#'            Y4 ~~ .301*Y5"
#' cfaHB(manmod,500,string=TRUE)
cfaHB <- function(model,n=NULL,plot=FALSE,string=FALSE){

  #If string, expect string (a la Shiny app)
  if(string){
    model=model
    n=n
  }else{
    #Use these functions to convert to string (input is a lavaan object)
    #Probably what we should expect for people using R
    #need 'n' first because otherwise model will overwrite
    n <- cfa_n(model)
    model <- cfa_lavmod(model)
  }

  if (unstandardized(model)>0){
    stop("dynamic Error: Your model has loadings greater than or equal to 1 (an impossible value). Please use standardized loadings.")
  }

  if (number_factor(model)<2){
    stop("dynamic Error: You entered a one-factor model.  Use cfaOne instead.")
  }

  if (defre(model,n)==0){
    stop("dynamic Error: It is impossible to add misspecifications to a just identified model.")
  }

  if ( nrow(multi_num_HB(model)) < (number_factor(model)-1)){
    stop("dynamic Error: There are not enough free items to produce all misspecification levels.")
  }

  #Create list to store outputs (table and plot)
  res <- list(input=as.list(environment),
              output=list())

  #Run simulation
  results <- multi_df_HB(model,n)

  #For each list element (misspecification) compute the cutoffs
  misspec_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_M=quantile(SRMR_M, c(.05,.1)),
                                                      RMSEA_M=quantile(RMSEA_M, c(.05,.1)),
                                                      CFI_M=quantile(CFI_M, c(.95,.9))))

  #For the true model, compute the cutoffs (these will all be the same - just need in list form)
  true_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_T=quantile(SRMR_T, c(.95,.9)),
                                                   RMSEA_T=quantile(RMSEA_T, c(.95,.9)),
                                                   CFI_T=quantile(CFI_T, c(.05,.1))))

  #Bind each of the misspecified cutoffs to the true cutoffs, listwise
  Table <- purrr::map(misspec_sum,~cbind(.,true_sum[[1]]) %>%
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
                           dplyr::mutate_at(c("SRMR_1","RMSEA_1","CFI_1"),list(dplyr::lead)) %>%
                           dplyr::slice(1) %>%
                           dplyr::mutate(SRMR=ifelse(is.character(SRMR),SRMR_1,"--"),
                                         RMSEA=ifelse(is.character(RMSEA),RMSEA_1,"--"),
                                         CFI=ifelse(is.character(CFI),CFI_1,"--"),
                                         SRMR=stringr::str_replace_all(as.character(SRMR),"0\\.","."),
                                         RMSEA=stringr::str_replace_all(as.character(RMSEA),"0\\.","."),
                                         CFI=stringr::str_replace_all(as.character(CFI),"0\\.",".")) %>%
                           dplyr::select(SRMR,RMSEA,CFI))

  #Still cleaning
  #Unlist Table
  Table_C <- purrr::map_dfr(Table,~dplyr::mutate(.,SRMR=stringr::str_replace_all(as.character(SRMR),"0\\.","."),
                                                 RMSEA=stringr::str_replace_all(as.character(RMSEA),"0\\.","."),
                                                 CFI=stringr::str_replace_all(as.character(CFI),"0\\.",".")))

  #Cleaning
  Table_C[seq(2,nrow(Table_C),by=2),] <- Row2

  #For row names
  num_fact <- (number_factor(model)-1)

  #Create row names for level
  Table_C$levelnum <- paste("Level", rep(1:num_fact,each=2))

  #Create row names for proportions
  Table_C$cut <- rep(c("95/5","90/10"))

  #Add cross-loading magnitude
  suppressMessages(mag <- multi_add_HB(model) %>%
                     tidyr::separate(V1,into=c("a","b","Magnitude","d","e"),sep=" ") %>%
                     select(Magnitude) %>%
                     mutate(Magnitude=as.numeric(Magnitude),
                            Magnitude=round(Magnitude,digits=3)) %>%
                     slice(rep(1:n,each=2)))

  #Clean cross-loading magnitude
  even <- seq_len(nrow(mag))%%2
  mag2 <- cbind(mag,even) %>%
    mutate(Magnitude=ifelse(even==0," ",Magnitude)) %>%
    mutate(Magnitude=stringr::str_replace_all(as.character(Magnitude),"0\\.",".")) %>%
    select(Magnitude)

  #Add to table
  Table_C <- cbind(Table_C,mag2)

  #Add rownames to final table
  Final_Table <- Table_C %>%
    tidyr::unite(Cut,levelnum,cut,sep=": ") %>%
    tibble::column_to_rownames(var='Cut')

  #Put into list
  res$output$Cutoffs <- Final_Table

  #If user selects plot = T
  if(plot){
    #For each list element (misspecification) compute the cutoffs
    misspec_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_M=quantile(SRMR_M, c(.05,.1)),
                                                        RMSEA_M=quantile(RMSEA_M, c(.05,.1)),
                                                        CFI_M=quantile(CFI_M, c(.95,.9))))

    #For the true model, compute the cutoffs (these will all be the same - just need in list form)
    true_sum <- purrr::map(results,~dplyr::summarise(.,SRMR_T=quantile(SRMR_T, c(.95,.9)),
                                                     RMSEA_T=quantile(RMSEA_T, c(.95,.9)),
                                                     CFI_T=quantile(CFI_T, c(.05,.1))))

    #Select just those variables and rename columns to be the same
    Misspec_dat <- purrr::map(results,~dplyr::select(.,SRMR_M:Type_M) %>%
                                `colnames<-`(c("SRMR","RMSEA","CFI","Model")))

    #Select just those variables and rename columns to be the same
    True_dat <- purrr::map(results,~dplyr::select(.,SRMR_T:Type_T) %>%
                             `colnames<-`(c("SRMR","RMSEA","CFI","Model")))

    #For each element in the list, bind the misspecified cutoffs to the true cutoffs
    #rbind doesn't work well with lists (needs do.call statement)
    plot <- lapply(seq(length(Misspec_dat)),function(x) dplyr::bind_rows(Misspec_dat[x],True_dat[x]))

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
    plots_combo <- lapply(seq(length(plot)),function(x) c(SRMR_plot[x],RMSEA_plot[x],CFI_plot[x]))

    #Add a collective legend and title with the level indicator
    plots <- lapply(seq(length(plots_combo)), function(x) patchwork::wrap_plots(plots_combo[[x]])+
                      plot_layout(guides = "collect")+
                      plot_annotation(title=paste("Level", x))
                    & theme(legend.position = 'bottom'))

    #Put into list
    res$output$Plots <- plots

  }

  #Create object (necessary for subsequent print statement)
  class(res) <- 'cfaHB'

  return(res)

}

#' @method print cfaHB
#' @param x cfaHB object
#' @param ... other print parameters
#' @rdname cfaHB
#' @export

#Print suppression/organization statement for list
#Needs same name as class, not function name
#Need to add ... param or will get error message in CMD check
print.cfaHB <- function(x,...){

  base::cat("Your DFI cutoffs: \n")
  base::print(x$output$Cutoffs)

  if(!is.null(x$output$Plots)){

    base::cat("\n The distributions for each level are in the Plots tab \n")
    base::print(x$output$Plots)
  }

  #Hides this function
  base::invisible()
}
