#' @title Deprecated: Fit index cutoffs for single level confirmatory factor analysis models
#'
#' @description Deprecated. If you want to replicate the cfaFit cutoffs for multi-factor models, use the function \code{\link{cfaHB}}.
#' The first level will be identical to the results from cfaFit. The cutoffs for one-factor models from cfaFit cannot be
#' replicated in other functions. Use \code{\link{cfaOne}} for one-factor CFA models.
#'
#' @importFrom tibble column_to_rownames
#' @import patchwork
#'
#' @param model Model statement written in \code{\link{lavaan}} \code{\link{model.syntax}} with standardized loadings.
#'     If using \code{\link{lavaan}}, you can get standardized loadings from \code{\link{standardizedSolution}}.
#' @param n Sample size (number).
#'
#' @return Fit index cutoffs for the SRMR, RMSEA, and CFI.
#' @export
#'
#' @examples
#' mod <- "F1 =~ .602*Y1 + .805*Y2 + .516*Y3 + .857*Y4
#'         F2 =~ .413*Y5 + -.631*Y6
#'         F1 ~~ .443*F2
#'         Y4 ~~ .301*Y5"
#' cfaFit(mod,500)
cfaFit <- function(model,n){

  if (unstandardized(model)>0){
    stop("dynamic Error: Your model has loadings greater than or equal to 1 (an impossible value). Please use standardized loadings.")
  }

  .Deprecated(msg="cfaFit is deprecated and will be defunct in the next version of dynamic.

  If you want to replicate the cfaFit cutoffs for multi-factor models, use the function cfaHB.  The first level will be identical to the results from cfaFit. The cutoffs for one-factor models from cfaFit cannot be replicated in other functions. Use cfaOne for one-factor CFA models.")

  results <- dynamic_fit(model,n)

  misspec_sum <- results %>%
    dplyr::summarise(SRMR_M=quantile(SRMR_M, c(.05,.1)),
                     RMSEA_M=quantile(RMSEA_M, c(.05,.1)),
                     CFI_M=quantile(CFI_M, c(.95,.9)))

  true_sum <- results %>%
    dplyr::summarise(SRMR_T=quantile(SRMR_T, c(.95,.9)),
                     RMSEA_T=quantile(RMSEA_T, c(.95,.9)),
                     CFI_T=quantile(CFI_T, c(.05,.1)))

  Table <- cbind(misspec_sum,true_sum) %>%
    dplyr::mutate(SRMR_R=base::round(SRMR_M,3),
                  RMSEA_R=base::round(RMSEA_M,3),
                  CFI_R=base::round(CFI_M,3),
                  SRMR=base::ifelse(SRMR_T<SRMR_M,SRMR_R,"NONE"),
                  RMSEA=base::ifelse(RMSEA_T<RMSEA_M,RMSEA_R,"NONE"),
                  CFI=base::ifelse(CFI_T>CFI_M,CFI_R,"NONE")) %>%
    dplyr::select(SRMR,RMSEA,CFI)

  T_R1 <- Table %>%
    dplyr::slice(1) %>%
    dplyr::rename(SRMR_1=SRMR,
                  RMSEA_1=RMSEA,
                  CFI_1=CFI)

  T_R2 <- Table%>%
    dplyr::slice(2) %>%
    dplyr::rename(SRMR_2=SRMR,
                  RMSEA_2=RMSEA,
                  CFI_2=CFI)

  T_C <- cbind(T_R1,T_R2) %>%
    dplyr::mutate(SRMR=base::ifelse(base::is.character(SRMR_1),SRMR_2,"--"),
                  RMSEA=base::ifelse(base::is.character(RMSEA_1),RMSEA_2,"--"),
                  CFI=base::ifelse(base::is.character(CFI_1),CFI_2,"--")) %>%
    dplyr::select(SRMR,RMSEA,CFI)

  Table_Final <- Table %>%
    dplyr::slice(1) %>%
    base::rbind(T_C) %>%
    dplyr::mutate(Cut=c("95/5","90/10")) %>%
    dplyr::select(Cut,SRMR,RMSEA,CFI) %>%
    tibble::column_to_rownames(var="Cut")

  Misspec_dat <- results %>%
    dplyr::select(SRMR_M:Type_M) %>%
    `colnames<-`(c("SRMR","RMSEA","CFI","Model"))

  True_dat <- results %>%
    dplyr::select(SRMR_T:Type_T) %>%
    `colnames<-`(c("SRMR","RMSEA","CFI","Model"))

  misspec_sum <- results %>%
    dplyr::summarise(SRMR_M=quantile(SRMR_M, c(.05,.1)),
                     RMSEA_M=quantile(RMSEA_M, c(.05,.1)),
                     CFI_M=quantile(CFI_M, c(.95,.9)))

  true_sum <- results %>%
    dplyr::summarise(SRMR_T=quantile(SRMR_T, c(.95,.9)),
                     RMSEA_T=quantile(RMSEA_T, c(.95,.9)),
                     CFI_T=quantile(CFI_T, c(.05,.1)))

  plot <- base::rbind(Misspec_dat,True_dat)

  SRMR_plot <- plot %>%
    ggplot(aes(x=SRMR,fill=Model))+
    geom_histogram(position="identity",
                   bins = 30,
                   alpha=.5)+
    scale_fill_manual(values=c("#E9798C","#66C2F5"))+
    geom_vline(aes(xintercept=misspec_sum$SRMR_M[1],
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
          legend.box = "vertical")

  RMSEA_plot <- plot %>%
    ggplot(aes(x=RMSEA,fill=Model))+
    geom_histogram(position="identity",
                   bins = 30,
                   alpha=.5)+
    scale_fill_manual(values=c("#E9798C","#66C2F5"))+
    geom_vline(aes(xintercept=misspec_sum$RMSEA_M[1],
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
          legend.box = "vertical")

  CFI_plot <- plot %>%
    ggplot(aes(x=CFI,fill=Model))+
    geom_histogram(position="identity",
                   bins = 30,
                   alpha=.5)+
    scale_fill_manual(values=c("#E9798C","#66C2F5"))+
    geom_vline(aes(xintercept=misspec_sum$CFI_M[1],
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
          legend.box = "vertical")

  #Plot_Final <- base::suppressMessages(lemon::grid_arrange_shared_legend(SRMR_plot, RMSEA_plot, CFI_plot,
                                                                         #ncol=3,nrow=1))
  Plot_Final <- SRMR_plot + RMSEA_plot + CFI_plot +
    patchwork::plot_layout(guides = "collect") & theme(legend.position = 'bottom')

  return(list(Table_Final, Plot_Final))
}

