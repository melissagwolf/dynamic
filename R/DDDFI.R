#' @title Direct Discrepancy Dynamic fit index (3DFI) cutoffs for arbitrary covariance structure models
#'
#' @description This function generates DFI cutoffs for any single group covariance structure model with
#' no mean structure or a saturated mean structure. It supports (a) any estimator supported by lavaan (e.g.
#' ML, MLR, WLSMV, ULSMV), (b) missing data, and (c) multiple data formats (normal, non-normal continuous,
#' categorical). The default argument is a singular argument: a \code{\link{lavaan}} object.The function
#' can also accommodate manual entry of the model statement and sample size. Some features require an original
#' dataset to be provided as well (e.g., missing data, categorical data). The app-based version of
#' this function can be found at \href{https://dynamicfit.app/}{dynamicfit.app}.
#'
#' @param model This can either be a \code{\link{lavaan}} object, OR a model statement written in \code{\link{lavaan}} \code{\link{model.syntax}}
#' with standardized estimates
#' @param scale Determines how data are simulated in the function. Options are "normal", "nonnormal", or "categorical". "normal" assumes multivariate
#' normality across all variables. "nonnormal" recreates distributions in an empirical dataset (to be provided by the user), assuming variables are continuous.
#' "categorical" simulates discrete data with the same proportions as an empirical dataset (to be provided by the user). With "categorical", mixed formats
#' are also supported and any variable with more than 9 categories is simulated from a normal distribution.
#' @param data The original data to which the model was applied. Not required if scale="normal". Otherwise, data is required.
#' @param manual If you entered a \code{\link{lavaan}} object, keep this set to FALSE.
#' If you manually entered standardized estimates and sample size, set this to TRUE.
#' @param reps The number of replications used in your simulation. This is set to 250 by default
#' @param n If you entered a \code{\link{lavaan}} object for model, leave this blank.
#' Otherwise, enter your sample size (numeric).
#' @param estimator Which estimator to use within the simulations (enter in quotes). The default is "ML".
#' @param MAD Mean Absolute Discrepancies to test in the simulation. Default is c(.038, .05, .06) to recreate traditional "Close", "Fair", "Mediocre" benchmarks
#' @param plot.dfi Displays simulated distributions of fit indices used to derive cutoffs for each level of misspecification.
#' @param plot.dist Displays distributions of simulated data (and empirical data, if provided) to assess fidelity of simulated data to empirical data
#' @param plot.discrepancy Displays distributions of simulated MAD values
#'
#' @import dplyr lavaan simstandard ggplot2 stringr MASS
#' @importFrom purrr map map_dfr map2
#' @importFrom tidyr unite extract
#' @importFrom patchwork plot_layout plot_annotation wrap_plots
#' @importFrom semTools bsBootMiss

#'
#' @author Daniel McNeish & Melissa G Wolf
#'
#' Maintainer: Daniel McNeish <dmcneish@asu.edu>
#'
#' @rdname DDDFI
#'
#' @return Direct Discrepancy Dynamic fit index (DFI) cutoffs for CFI, RMSEA, and RMSEA 90% CI.
#' @export

###########################################################
# Bookkeeping functions
# Used to format input (from lavaan object or manual entry)
# Or to catch input errors/format output
###########################################################

#Function to strip estimates from model statement
#Used to apply fitted model to simulated data
cleanmodel <- function(model){

  suppressMessages(  model %>%
                       lavaan::lavaanify(fixed.x = FALSE) %>%
                       dplyr::filter(lhs != rhs) %>%
                       dplyr::filter(op != "~1") %>%
                       dplyr::filter(op != "|") %>%
                       dplyr::group_by(lhs, op) %>%
                       mutate(rhs=case_when(ustart==0 ~paste0(ustart,"*",rhs),
                                            ustart!=0 ~ rhs)) %>%
                       dplyr::summarise(rhs = paste(rhs, collapse = " + ")) %>%
                       dplyr::arrange(dplyr::desc(op)) %>%
                       tidyr::unite("l", lhs, op, rhs, sep = " ") %>%
                       dplyr::pull(l))

}

#function checking for inadmissible estimates
#used in warning messages
unstandardized <- function(model){

  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(lhs != rhs)

  one_plus <- lav_file %>%
    dplyr::filter(ustart >= 1 | ustart <= -1) %>%
    base::nrow()

  return(one_plus)
}

#note-- do we need a function for this if it's one line?
#function to extract sample size from lavaan object
cfa_n <- function(model){

  #Extract n from lavaan object
  #Warning message to hide warning when someone enters a non-lavaan object (error message will display instead)
  #Only need it here because this is the first argument in cfaHB and cfaOne
  n <- base::unlist(model@SampleStats@nobs)
  return(n)
}

#extract model statement from lavaan object
cfa_lavmod <- function(model){

  #Extract standardized solution from lavaan object
  lav <- lavaan::standardizedSolution(model)

  #Create model statement
  ss_mod <- suppressMessages(lav %>%
                               dplyr::filter(lhs != rhs) %>%
                               dplyr::group_by(lhs,op) %>%
                               dplyr::filter(op != "~1") %>%
                               dplyr::filter(op != "|") %>%
                               dplyr::select(lhs,op,rhs,est.std) %>%
                               dplyr::mutate(est.std=round(est.std,digits=4)) %>%
                               dplyr::summarise(rhs=paste(est.std,"*",rhs,collapse=" + ")) %>%
                               dplyr::arrange(desc(op)) %>%
                               tidyr::unite("mod",lhs,op,rhs,sep="") %>%
                               dplyr::pull(mod))

  #Collapse into one string because my other functions expect that
  mod <- base::paste(ss_mod, sep="", collapse="\n")

  return(mod)
}


# function to remove thresholds from lavaan-input model statement (if categorical =T)
# thresholds not needed to generate data, but are included with output
modelWithNum <- function(model){

  suppressMessages(model %>%
                     lavaan::lavaanify(fixed.x=FALSE) %>%
                     dplyr::filter(lhs != rhs) %>%
                     dplyr::filter(op != "~1") %>%
                     dplyr::filter(op != "|") %>%
                     dplyr::group_by(lhs,op) %>%
                     tidyr::unite("xx",ustart,rhs, sep="*")  %>%
                     dplyr::summarise(rhs = paste(xx, collapse = " + "))%>%
                     dplyr::arrange(dplyr::desc(op))%>%
                     tidyr::unite("l", lhs, op, rhs, sep = " ") %>%
                     summarise(l=paste(l, collapse="\n")) %>%
                     dplyr::pull(l)
  )
}

#function to collate fit index values across replications
#part of optional output
fit_data <- function(df_results){

  #Create beginning of variable name for each
  dat_name <- base::rep(c("SRMR_L","RMSEA_L","CFI_L","Type_L"),2)

  #Create vector of 0's for the True model
  dat_0 <- base::rep(0,4)

  #Get number of levels of misspecification
  dat_lev <- base::length(df_results)

  #Create combo of Level #'s and 0's to merge with variable names (in list form)
  dat_num <- base::list()
  for (i in 1:dat_lev){
    output <- c(base::rep(i,4),dat_0)
    dat_num[[i]] <- output
  }

  #Combine variable name with level # (in list form)
  var_names <- base::lapply(dat_num, function(x){
    base::paste0(dat_name,x)
  })

  #Rename variables in dataset list
  dat_revised <- base::lapply(base::seq_along(df_results), function(x){
    base::colnames(df_results[[x]]) <- var_names[[x]]
    #not sure why I need to mention it again but I do
    df_results[[x]]
  })

  #Combine into one dataset
  df_renamed <- do.call(base::cbind.data.frame,dat_revised)

  #Remove the duplicate L0 info
  df_renamed2 <- df_renamed[!base::duplicated(base::colnames(df_renamed))]

  #Remove the "type" variable (unnecessary)
  df_renamed3 <- df_renamed2[, -base::grep("Type",base::colnames(df_renamed2))]

  #Reorder variables
  df_renamed4 <- df_renamed3 %>%
    dplyr::relocate(base::order(base::colnames(df_renamed3))) %>%            #gets numbers in order
    dplyr::relocate(dplyr::starts_with("SRMR"),dplyr::starts_with("RMSEA"))    #gets fit indices in order

  return(df_renamed4)
}

###############################################################################
# Specific functions from other packages, lightly modified to serve our purpose
# Also reduces dependencies and susceptibility to changes in other packages
###############################################################################

#There is no build in function for the LKJ distribution
#Manually write out LKJ distribution (taken from other R package, XXX)
rlkjcorr <- function (n, K, eta = 1) {
  stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
  stopifnot(eta > 0)
  #if (K == 1) return(matrix(1, 1, 1))

  f <- function() {
    alpha <- eta + (K - 2)/2
    r12 <- 2 * rbeta(1, alpha, alpha) - 1
    R <- matrix(0, K, K) # upper triangular Cholesky factor until return()
    R[1,1] <- 1
    R[1,2] <- r12
    R[2,2] <- sqrt(1 - r12^2)
    if(K > 2) for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m / 2, alpha)

      # Draw uniformally on a hypersphere
      z <- rnorm(m, 0, 1)
      z <- z / sqrt(crossprod(z)[1])

      R[1:m,m+1] <- sqrt(y) * z
      R[m+1,m+1] <- sqrt(1 - y)
    }
    return(crossprod(R))
  }
  R <- replicate( n , f() )
  if ( dim(R)[3]==1 ) {
    R <- R[,,1]
  } else {
    # need to move 3rd dimension to front, so conforms to array structure that Stan uses
    R <- aperm(R,c(3,1,2))
  }
  return(R)
}

##function to cut data based on quantiles
##used for generated discrete data treated as continuous
#taken from gtools version 3.9.4
quantcut <- function(x, q = 4, na.rm = TRUE, ...) {
  if (length(q) == 1) {
    q <- seq(0, 1, length.out = q + 1)
  }

  quant <- quantile(x, q, na.rm = na.rm)
  dups <- duplicated(quant)
  if (any(dups)) {
    flag <- x %in% unique(quant[dups])
    retval <- ifelse(flag,
                     paste("[",
                           as.character(x),
                           "]",
                           sep = ""
                     ),
                     NA
    )
    uniqs <- unique(quant)

    # move cut points over a bit...
    reposition <- function(cut) {
      flag <- x >= cut
      if (sum(flag, na.rm = na.rm) == 0) {
        return(cut)
      } else {
        return(min(x[flag], na.rm = na.rm))
      }
    }

    newquant <- sapply(uniqs, reposition)
    retval[!flag] <- as.character(cut(x[!flag],
                                      breaks = newquant,
                                      include.lowest = TRUE, ...
    ))

    levs <- unique(retval[order(x)]) # ensure factor levels are
    # properly ordered
    retval <- factor(retval, levels = levs)

    ## determine open/closed interval ends
    mkpairs <- function(x) { # make table of lower, upper
      sapply(
        x,
        function(y) if (length(y) == 2) y[c(2, 2)] else y[2:3]
      )
    }
    pairs <- mkpairs(strsplit(levs, "[^0-9+\\.\\-]+"))
    rownames(pairs) <- c("lower.bound", "upper.bound")
    colnames(pairs) <- levs

    closed.lower <- rep(F, ncol(pairs)) # default lower is open
    closed.upper <- rep(T, ncol(pairs)) # default upper is closed
    closed.lower[1] <- TRUE # lowest interval is always closed

    for (i in 2:ncol(pairs)) { # open lower interval if above singlet
      if (pairs[1, i] == pairs[1, i - 1] && pairs[1, i] == pairs[2, i - 1]) {
        closed.lower[i] <- FALSE
      }
    }

    for (i in 1:(ncol(pairs) - 1)) { # open upper inteval if below singlet
      if (pairs[2, i] == pairs[1, i + 1] && pairs[2, i] == pairs[2, i + 1]) {
        closed.upper[i] <- FALSE
      }
    }

    levs <- ifelse(pairs[1, ] == pairs[2, ],
                   pairs[1, ],
                   paste(ifelse(closed.lower, "[", "("),
                         pairs[1, ],
                         ",",
                         pairs[2, ],
                         ifelse(closed.upper, "]", ")"),
                         sep = ""
                   )
    )
    levels(retval) <- levs
  }
  else {
    retval <- cut(x, quant, include.lowest = TRUE, ...)
  }
  return(retval)
}

#############################################################
# Functions to Generate Data and Fit Models to Simulated Data
#############################################################

#formerly multi_fit
#Created discrepancy matrices
#generated misspecified data
#fit original model to misspecified data
#collated fit indices into a data frame
miss_fit <- function(model,data,n,reps,estimator,MAD,scale){

  #strip estimates from model statement
  mod <- cleanmodel(model)

  # if categorical is True, remove thresholds from model statement as well
  if((scale %in% c("categorical"))) {
    model<-modelWithNum(model)
  }

  #count the number of factors in the model (needed to get dimension of matrix only for observed variables)
  factors <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(lhs != rhs)%>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #simulate model-implied (polychoric) correlation matrix for observed variables
  dat <- simstandard::sim_standardized_matrices(model)
  x<-dat$Correlations$R
  l<-nrow(x)-nrow(factors)
  a<-x[1:l,1:l]
  a<-a[order(rownames(a)),order(colnames(a))]

  r<-reps
  n <- base::min(n,5000)

  #create lists to house simulated, misspecified data

  all_data_misspec<-vector(mode="list",length=length(MAD))
  all_Q<-vector(mode="list",length=length(MAD))

  #if scale="normal", simulate everything from multivariate normal, no need for original data
  if (scale %in% c("normal")){

    set.seed(97)
    mu=rep(0,l)
    all_data_misspec<-vector(mode="list",length=length(MAD))
    all_Q<-vector(mode="list",length=length(MAD))

    #create discrepancy matrices and create misspecified data
    for (m in 1:length(MAD)) {
      M<-vector(mode="list",length=r)
      L<-vector(mode="list",length=r)
      D<-vector(mode="list",length=r)
      Q<-vector(mode="list",length=r)

      #select MAD from list
      mad<-MAD[m]
      #find row from LKJ grid that most closely match desired MAD
      row<-lkj[lkj$dim==l,]
      #select eta value to use for LKJ distribution
      eta<-as.numeric(row[which.min(abs(row$mean - mad)),]$eta)

      for(i in 1:r)
      {

        #root of discrepancy matrix plus root of model-implied matrix (avoids possible non-positive definite issues)
        M[[i]] <- chol(a)+chol(rlkjcorr(1,l,eta=eta))-diag(l)
        #full data generation correlation matrix (includes misspecifications)
        L[[i]] <- cov2cor(t(M[[i]])%*%M[[i]])
        #simulate data from full matrix
        D[[i]]<-as.data.frame(mvrnorm(n=n, mu=mu, Sigma=L[[i]]))
        #record replication number
        D[[i]]$rep <-i
        #track discrepancy of each element (optional part of output)
        Q[[i]]<-a-L[[i]]
      }
      #rename objects with suffix for misspecification level
      all_data_misspec[[m]]<-D
      assign(paste0("m",m),M)
      assign(paste0("l",m),L)
      assign(paste0("data_mis",m),D)

      #save all individual discrepancies (for optional output)
      all_Q[[m]]<-Q
    }
  }

  #if scale="nonnormal", used Fleishman method using original data to get skew and kurtosis
  if (scale %in% c("nonnormal")){

    unique<-lengths(lapply(data[,colnames(a)], unique))
    #flag any variable with between 2 and 7 categories are categorical/Likert
    probLik <- (1< unique & unique <10)
    probLik1 <- t(as.data.frame(unique[probLik]))
    #save names of likely categorical/likert variables (to be transformed later)
    likertnames<-colnames(probLik1)

    #create empty matrix for proportions in each category,g
    g<-list()
    #data only with discrete items
    data1<-data[,likertnames]

    #rescale so that minimum value is always ==1
    #needed to properly index computations below
    #will be scaled back at the end after loop indexes are not longer needed
    d2<-matrix(sapply(data1, function(x) min(x, na.rm=T)-1), nrow=1, ncol=ncol(data1))
    d3<-matrix(rep(d2,each=nrow(data1)), nrow=nrow(data1), ncol=ncol(data1))
    d3l<-matrix(rep(d2,each=(r*nrow(data1))), nrow=(r*nrow(data1)), ncol=ncol(data1))
    colnames(d3)<-colnames(data1)
    colnames(d3l)<-colnames(data1)
    data1<-data1-d3

    #loop through all discrete variables in fitted model
    for (i in 1:length(likertnames)){
      #setup list element of all 0s, to be replaced
      g1<-rep(0,max(data1[,i], na.rm=T)-1)
      g[[i]]<-g1
      #loop through all categories for each varaible
      for (h in (min(data1[,i], na.rm=T)-1):(max(data1[,i], na.rm=T))){
        #proportion of responses at or below each category, only count non-missing in the denominator
        g[[i]][h+1]<-sum(table(data1[,i])[names(table(data1[,i]))<=h])/sum(!is.na(data1[,i]))
      }
    }


    #calculate EM mean/cov?
    dummymod<-sem(mod,meanstructure=T,data=data[,colnames(a)],missing="ML",do.fit=F)
    #em<-lavInspect(dummymod,what="sampstat")
    em<-lavInspect(dummymod,what="sampstat")
    mu<-em$mean[colnames(a)]
    x<-chol(em$cov[colnames(a),colnames(a)])
    xx<-t(x)%*%x
    xxx<-xx[colnames(a),colnames(a)]
    S.inv.sqrt <- solve(chol(xxx))
    ##data for transforming
    ddd<-as.matrix(data[,colnames(a)])
    d<-ddd
    #missing data patterns
    na<-is.na(data[,colnames(a)])
    #groups of missing data
    nagrp<-unique(na)

    for (m in 1:length(MAD)) {
      M<-vector(mode="list",length=r)
      L<-vector(mode="list",length=r)
      D<-vector(mode="list",length=r)
      Q<-vector(mode="list",length=r)
      #select MAD from list
      mad<-MAD[m]
      #find row from LKJ grid that most closely match desired MAD
      row<-lkj[lkj$dim==l,]
      #select eta value to use for LKJ distribution
      eta<-as.numeric(row[which.min(abs(row$mean - mad)),]$eta)

      for(i in 1:r){

        ll<-list()
        gg<-vector()

        #root of discrepancy matrix plus root of model-implied matrix (avoids possible non-positive definite issues)
        M[[i]] <- chol(a)+chol(rlkjcorr(1,l,eta=eta))-diag(l)
        #full data generation correlation matrix (includes misspecifications)
        L[[i]] <- cov2cor(t(M[[i]])%*%M[[i]])
        #BS-bootstrap transformation
        sigma.sqrt <- chol(L[[i]])

        for (e in 1:nrow(nagrp)){
          for (f in 1:nrow(data)) {
            gg[f] <- all(na[f,]==nagrp[e,])
          }
          ll[[e]]<-gg
          dataz<-matrix(ddd[ll[[e]],nagrp[e,]==F],nrow=sum(ll[[e]]==T),ncol=length(names(nagrp[e,nagrp[e,]==F])))
          colnames(dataz)<-names(nagrp[e,nagrp[e,]==F])

          zz<-matrix(1,nrow(dataz),1)%*%matrix(1,1,ncol(dataz))+(dataz-matrix(1,nrow(dataz),1)%*%t(mu[colnames(dataz)]))%*%S.inv.sqrt[colnames(dataz),colnames(dataz)]%*%sigma.sqrt[colnames(dataz),colnames(dataz)]

          d[ll[[e]],nagrp[e,]==F]<-zz
        }

        #discretize discrete items
        #for(k in 1:length(likertnames)){
        #  d[,likertnames[k]]<-as.integer(quantcut(d[,likertnames[k]],q=g[[k]]))
        #}

        #add d3 back to put categorical items onto original scale
        d[,likertnames]<- d[,likertnames]+d3[,likertnames]
        #save data
        D[[i]]<-as.data.frame(d)
        D[[i]]$rep<-i

        #track discrepancy of each element (optional part of output)
        Q[[i]]<-a-L[[i]]
      }

      #rename objects with suffix for misspecification level
      all_data_misspec[[m]]<-D
      assign(paste0("m",m),M)
      assign(paste0("l",m),L)
      assign(paste0("data_mis",m),D)

      #save all individual discrepancies (for optional output)
      all_Q[[m]]<-Q
    }
  }

  #median(abs(L[[7]]-cor(D[[7]][,colnames(a)],use="complete.obs")))
  #mean(abs(a[colnames(a),colnames(a)]-cor(D[[8]][,colnames(a)],use="complete.obs")))
  #mean(abs(a[colnames(a),colnames(a)]-L[[10]]))

  #if scale="categorical", use data to figure out which variables are discrete and what proportions should be
  if (scale %in% c("categorical")){

    mu=rep(0,l)
    set.seed(97)
    all_data_misspec<-vector(mode="list",length=length(MAD))
    all_Q<-vector(mode="list",length=length(MAD))

    for (m in 1:length(MAD)) {
      M<-vector(mode="list",length=r)
      L<-vector(mode="list",length=r)
      D<-vector(mode="list",length=r)
      Q<-vector(mode="list",length=r)

      #select MAD from list
      mad<-MAD[m]
      #find row from LKJ grid that most closely match desired MAD
      row<-lkj[lkj$dim==l,]
      #select eta value to use for LKJ distribution
      eta<-as.numeric(row[which.min(abs(row$mean - mad)),]$eta)

      for(i in 1:r){

        #root of discrepancy matrix plus root of model-implied matrix (avoids possible non-positive definite issues)
        M[[i]] <- chol(a)+chol(rlkjcorr(1,l,eta=eta))-diag(l)
        #full data generation correlation matrix (includes misspecifications)
        L[[i]] <- cov2cor(t(M[[i]])%*%M[[i]])
        #simulate data from full matrix
        D[[i]]<-as.data.frame(mvrnorm(n=n, mu=mu, Sigma=L[[i]]))
        #record replication number
        D[[i]]$rep <-i
        #track discrepancy of each element (optional part of output)
        Q[[i]]<-a-L[[i]]
      }

      #rename objects with suffix for misspecification level
      all_data_misspec[[m]]<-D
      assign(paste0("m",m),M)
      assign(paste0("l",m),L)
      assign(paste0("data_mis",m),D)

      #save all individual discrepancies (for optional output)
      all_Q[[m]]<-Q
    }

    #number of columns, including rep counter
    last<-ncol(as.data.frame(all_data_misspec[[1]][[1]]))
    #get names of all variables in the model
    names<-colnames(a)
    #count number of unique values per variable
    unique<-lengths(lapply(data[,names], unique))
    #flag any variable with between 2 and 7 categories are categorical/Likert
    probLik <- (1< unique & unique <10)
    probLik1 <- t(as.data.frame(unique[probLik]))
    #save names of likely categorical/likert variables (to be transformed later)
    likertnames<-colnames(probLik1)

    #number of continuous variables
    n_cont<-ncol(a)-length(likertnames)
    #names of continuous variables
    contnames<-setdiff(rownames(a),likertnames)

    #create empty matrix for proportions in each category,g
    g<-list()
    #data only with discrete items
    data1<-data[,likertnames]

    #rescale so that minimum value is always ==1
    #needed to properly index computations below
    #also needed to avoid 0s because there is multiplication involved
    #will be scaled back at the end after loop indexes are not longer needed
    d2<-matrix(sapply(data1, function(x) min(x, na.rm=T)-1), nrow=1, ncol=ncol(data1))
    d3<-matrix(rep(d2,each=nrow(data1)), nrow=nrow(data1), ncol=ncol(data1))
    colnames(d3)<-colnames(data1)
    data1<-data1-d3

    #create empty matrix for proportions in each category (g) and pseudo-threshold corresponding to that proportion (p)
    g<-matrix(nrow=(max(data1,na.rm=T)-min(data1, na.rm=T)), ncol=ncol(data1))
    p<-g

    #loop through all variables in fitted model
    for (i in 1:ncol(data1)){
      #loop through all categories for each variable
      for (h in min(data1[,i],na.rm=T):(max(data1[,i],na.rm=T)-1)){
        #proportion of responses at or below each category, only count non-missing in the denominator
        g[h,i]<-sum(table(data1[,i])[names(table(data1[,i]))<=h])/sum(!is.na(data1[,i]))
        #inverse standard normal to determine pseudo-thresholds
        p[h,i]<-qnorm(c(g[h,i]))
        dp<-as.data.frame(p)
      }
    }
    #Assign column names so it's clear which thresholds go with which variable
    colnames(dp)<-colnames(data1)
    a2<-colnames(dp)

    #loop through MAD values
    for (x in 1:length(MAD))
    {
      #loop through replications
      for (xx in 1:r)
      {
        #loop through pseduo-thresholds
        for (i in 1:ncol(dp))
        {# first loop transforms highest category (important for binary variables)
          u<- as.numeric(sum(!is.na(dp[,i])))
          all_data_misspec[[x]][[xx]]<-all_data_misspec[[x]][[xx]] %>%
            dplyr::mutate(!!a2[i] := case_when(
              !!rlang::sym(a2[i])  <= dp[1,i] ~100,
              !!rlang::sym(a2[i]) >   dp[u,i] ~100*(u+1),
              TRUE ~ !!rlang::sym(a2[i]))
            )
          #second loop transforms all lower categories (for any ordinal/Likert variables)
          if(sum(!is.na(dp[,i])) > 1){
            for (j in 1:sum(!is.na(dp[,i]))) {
              all_data_misspec[[x]][[xx]]<-all_data_misspec[[x]][[xx]] %>%
                dplyr::mutate(!!a2[i] := case_when(
                  between( !!rlang::sym(a2[i]) , dp[j,i], dp[j+1,i]) ~ as.numeric(100*(j+1)),
                  TRUE~ !!rlang::sym(a2[i]))
                )
            }
          }
        }

        #loops multiply by 100 to avoid overwriting MVN data with values above 1
        #divide by 100 to put things back onto original metric
        #add d3 back to put simulated data back onto original scale (lowest values does not have to be 1)
        all_data_misspec[[x]][[xx]][,likertnames] <-(all_data_misspec[[x]][[xx]][,likertnames]/100) + d3[,likertnames]

        #mimic missing data pattern from original data
        #first, remove replication counter from last column
        temp_na<-all_data_misspec[[x]][[xx]][,-last]
        #identify which cells of original data matrix are missing
        temp_na[is.na(data[names])]<-NA
        #replace simulated values with NA if original is missing
        all_data_misspec[[x]][[xx]]<-temp_na
        #restore replication counter
        all_data_misspec[[x]][[xx]]$rep<-xx
      }
    }
  }

  #use "estimator =" option to figure out which corrected/scaled fit index is needed
  # also skip calculating SEs to speed up simulations, if estimator allows for it

  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("rmsea.ci.upper.scaled","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("rmsea.ci.upper.robust","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("rmsea.ci.upper","rmsea","cfi")
  }

  if((estimator=="ML" | estimator=="MLR") & scale %in% "nonnormal"){
    missing<-"ML"
  } else {
    missing<-"listwise"
  }

  #fit model to all simulated datasets
  #if categorical = T, treats categorical variables in "likertnames" as ordered categorical variables
  if(scale %in% c("categorical")){
    misspec_cfa <- purrr::map(all_data_misspec, function(x) purrr::map(x, function(y) lavaan::cfa(model = mod, data=y, estimator=estimator,std.lv=TRUE,se=se, ordered = likertnames, check.gradient=FALSE,
                                                                                                  check.post=FALSE,
                                                                                                  check.vcov=FALSE,
                                                                                                  control=list(rel.tol=.001))))
  }

  #if categorical ==F, then just treat all variables as continuous
  if(!(scale %in% c("categorical"))){
    misspec_cfa <- purrr::map(all_data_misspec, function(x) purrr::map(x, function(y) lavaan::cfa(model = mod, data=y, estimator=estimator,std.lv=TRUE,se=se,check.gradient=FALSE,
                                                                                                  check.post=FALSE,
                                                                                                  check.vcov=FALSE,
                                                                                                  control=list(rel.tol=.001))))
  }

  #Extract fit stats from each rep into a data frame
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, ind)) %>%
                                  `colnames<-`(c("RMSEA_CI_UPPER_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))

  #reset the seed
  set.seed(NULL)

  #create list for different outcomes
  miss<-list()

  #object with fit index values
  miss$misspec_fit_sum<-misspec_fit_sum
  #object with discrepancies tested
  miss$all_Q<-all_Q
  #object with all simulated data
  miss$Data<-all_data_misspec

  return(miss)

}

#formerly true_fit_multi
#generate data consistent with the original model
#fit orginal model to correct/consistent data
#collated fit indices into a data frame

true_fit<- function(model,data,n,reps, estimator, MAD,scale){

  #strip estimates from model statement
  mod <- cleanmodel(model)

  # if categorical==T, remove thresholds from model statement as well
  if(scale %in% c("categorical")) {
    model<-modelWithNum(model)
  }
  #count the number of factors in the model (needed to get dimension of matrix only for observed variables)
  factors <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(lhs != rhs)%>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #simulate model-implied (polychoric) correlation matrix for observed variables
  dat <- simstandard::sim_standardized_matrices(model)
  x<-dat$Correlations$R
  l<-nrow(x)-nrow(factors)
  a<-x[1:l,1:l]
  a<-a[order(rownames(a)),order(colnames(a))]
  diag(a)=1

  r<-reps

  #Use max sample size of 5000
  n <- base::min(n,5000)

  #Set seed
  set.seed(649364)

  #if scale="normal", simulate everything from multivariate normal, no need for original data
  if (scale %in% c("normal")){

    #set mean vector to 0
    mu=rep(0,l)

    #create list to house simulated data
    data_t<-list()

    #simulate data that are consistent with the original model's implied correlation matrix
    for(i in 1:r)
    {
      data_t[[i]]<-as.data.frame(mvrnorm(n=n, mu=mu, Sigma=a))
      data_t[[i]]$rep<-i
    }

  }

  if (scale %in% c("nonnormal")){

    unique<-lengths(lapply(data[,colnames(a)], unique))
    #flag any variable with between 2 and 7 categories are categorical/Likert
    probLik <- (1< unique & unique <10)
    probLik1 <- t(as.data.frame(unique[probLik]))
    #save names of likely categorical/likert variables (to be transformed later)
    likertnames<-colnames(probLik1)

    #create empty matrix for proportions in each category,g
    g<-list()
    #data only with discrete items
    data1<-data[,likertnames]

    #rescale so that minimum value is always ==1
    #needed to properly index computations below
    #will be scaled back at the end after loop indexes are not longer needed
    d2<-matrix(sapply(data1, function(x) min(x, na.rm=T)-1), nrow=1, ncol=ncol(data1))
    d3<-matrix(rep(d2,each=nrow(data1)), nrow=nrow(data1), ncol=ncol(data1))
    d3l<-matrix(rep(d2,each=(r*nrow(data1))), nrow=(r*nrow(data1)), ncol=ncol(data1))
    colnames(d3)<-colnames(data1)
    colnames(d3l)<-colnames(data1)
    data1<-data1-d3

    #loop through all discrete variables in fitted model
    for (i in 1:length(likertnames)){
      #setup list element of all 0s, to be replaced
      g1<-rep(0,max(data1[,i], na.rm=T)-1)
      g[[i]]<-g1
      #loop through all categories for each varaible
      for (h in (min(data1[,i], na.rm=T)-1):(max(data1[,i], na.rm=T))){
        #proportion of responses at or below each category, only count non-missing in the denominator
        g[[i]][h+1]<-sum(table(data1[,i])[names(table(data1[,i]))<=h])/sum(!is.na(data1[,i]))
      }
    }

    dummymod<-sem(mod,meanstructure=T,data=data[,colnames(a)],missing="ML",do.fit=F)
    em<-lavInspect(dummymod,what="sampstat")
    mu<-em$mean[colnames(a)]
    data1<-data[,colnames(a)]

    data_t<-semTools::bsBootMiss(Sigma =a, Mu=mu, rawData = data1, nBoot=reps,bootSamplesOnly  = TRUE)
    #for (u in 1:r){
    #for(k in 1:length(likertnames)){
    #data_t[[u]][,likertnames[k]]<-as.integer(quantcut(data_t[[u]][,likertnames[k]],q=g[[k]]))
    #}
    #data_t[[u]][,likertnames]<- data_t[[u]][,likertnames]+d3[,likertnames]
    #data_t[[u]]<-data_t[[u]][,1:ncol(a)]
    #data_t[[u]]$rep<-u
    #}
  }

  if (scale %in% c("categorical")){
    #Set seed
    set.seed(649364)
    #set mean vector to 0
    mu=rep(0,l)
    #create list to house simulated data
    data_t<-list()

    #simulate data that are consistent with the original model's implied correlation matrix
    for(i in 1:r)
    {
      data_t[[i]]<-as.data.frame(mvrnorm(n=n, mu=mu, Sigma=a))
      data_t[[i]]$rep<-i
    }

    #number of columns, including rep counter
    last<-ncol(as.data.frame(data_t[[1]]))
    #get names of all variables in the model
    names<-colnames(a)
    #count number of unique values per variable
    unique<-lengths(lapply(data[,colnames(a)], unique))
    #flag any variable with between 2 and 7 categories are categorical/Likert
    probLik <- (1< unique & unique <10)
    probLik1 <- t(as.data.frame(unique[probLik]))
    #save names of likely categorical/likert variables (to be transformed later)
    likertnames<-colnames(probLik1)

    #number of continuous variables
    n_cont<-ncol(a)-length(likertnames)
    #names of continuous variables
    contnames<-setdiff(rownames(a),likertnames)

    #create empty matrix for proportions in each category,g
    g<-list()
    #data only with discrete items
    data1<-data[,likertnames]

    #rescale so that minimum value is always ==1
    #needed to properly index computations below
    #also needed to avoid 0s because there is multiplication involved
    #will be scaled back at the end after loop indexes are not longer needed
    d2<-matrix(sapply(data1, function(x) min(x, na.rm=T)-1), nrow=1, ncol=ncol(data1))
    d3<-matrix(rep(d2,each=nrow(data1)), nrow=nrow(data1), ncol=ncol(data1))
    colnames(d3)<-colnames(data1)
    data1<-data1-d3

    #create empty matrix for proportions in each category (g) and pseudo-threshold corresponding to that proportion (p)
    g<-matrix(nrow=(max(data1,na.rm=T)-min(data1, na.rm=T)), ncol=ncol(data1))
    p<-g

    #loop through all variables in fitted model
    for (i in 1:ncol(data1)){
      #loop through all categories for each varaible
      for (h in min(data1[,i], na.rm=T):(max(data1[,i], na.rm=T)-1)){
        #proportion of responses at or below each category, only count non-missing in the denominator
        g[h,i]<-sum(table(data1[,i])[names(table(data1[,i]))<=h])/sum(!is.na(data1[,i]))
        #inverse standard normal to determine pseudo-thresholds
        p[h,i]<-qnorm(c(g[h,i]))
        dp<-as.data.frame(p)
      }
    }
    #inverse standard normal to determine pseudo-thresholds
    colnames(dp)<-colnames(data1)
    a2<-colnames(dp)

    #loop through replications
    for (xx in 1:r)
    {
      #loop through pseduo-thresholds
      for (i in 1:ncol(dp))
      {# first loop transforms highest category (important distinction for binary variables)
        u<- as.numeric(sum(!is.na(dp[,i])))
        data_t[[xx]]<-data_t[[xx]] %>%
          dplyr::mutate(!!a2[i] := case_when(
            !!rlang::sym(a2[i])  <= dp[1,i] ~100,
            !!rlang::sym(a2[i]) >   dp[u,i] ~100*(u+1),
            TRUE ~ !!rlang::sym(a2[i]))
          )
        #second loop transforms all lower categories
        if(sum(!is.na(dp[,i])) > 1){
          for (j in 1:sum(!is.na(dp[,i]))) {
            data_t[[xx]]<-data_t[[xx]] %>%
              dplyr::mutate(!!a2[i] := case_when(
                between( !!rlang::sym(a2[i]) , dp[j,i], dp[j+1,i]) ~ as.numeric(100*(j+1)),
                TRUE~ !!rlang::sym(a2[i]))
              )
          }
        }
      }

      #loops multiply by 100 to avoid overwriting MVN data with values above 1
      #divide by 100 to put things back onto original metric
      #add d3 back to put simulated data back onto original scale (lowest values does not have to be 1)
      data_t[[xx]][,likertnames] <-(data_t[[xx]][,likertnames]/100) + d3[,likertnames]

      #######
      #mimic missing data pattern from original data
      ########

      #first, remove replication counter from last column
      #mimic missing data pattern from original data
      temp_na<-data_t[[xx]][,-last]
      #identify which cells of original data matrix are missing
      temp_na[is.na(data[names])]<-NA
      #replace simulated values with NA if original is missing
      data_t[[xx]]<-temp_na
      #restore replication counter
      data_t[[xx]]$rep<-xx
    }
  }

  # use "estimator=" option to determine which type of fit index to track (standard, normal, scaled)
  # also skip calculating SEs to speed up simulations, if estimator allows for it
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("rmsea.ci.upper.scaled","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("rmsea.ci.upper.robust","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("rmsea.ci.upper","rmsea","cfi")
  }

  if((estimator=="ML" | estimator=="MLR") & scale %in% "nonnormal"){
    missing<-"ML"
  } else {
    missing<-"listwise"
  }

  #fit model to all simulated datasets
  #if categorical = T, treats categorical variables in "likertnames" as ordered categorical variables
  if(scale %in% c("categorical")){
    true_cfa <- purrr::map(data_t, function(x) lavaan::cfa(model = mod, data=x, std.lv=TRUE, ordered=likertnames,estimator=estimator,se=se, check.gradient=FALSE,
                                                           check.post=FALSE,
                                                           check.vcov=FALSE,
                                                           control=list(rel.tol=.001)))
  }

  #if categorical ==F, then just treat all variables as continuous
  if(!(scale %in% c("categorical"))){
    true_cfa <- purrr::map(data_t, function(x) lavaan::cfa(model = mod, data=x, std.lv=TRUE,estimator=estimator,se=se,check.gradient=FALSE,
                                                           check.post=FALSE,
                                                           check.vcov=FALSE,
                                                           control=list(rel.tol=.001)))
  }


  #Extract fit stats from each rep into a data frame
  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., ind)) %>%
    `colnames<-`(c("RMSEA_CI_UPPER_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  #create list for different outcomes
  true<-list()

  #object with fit index values
  true$true_fit_sum<-true_fit_sum
  #objects with discrepancies tested
  true$all_D<-data_t

  return(true)
}

#formerly multi_df
#Function to combine all indices into one dataframe
#Will used to created distributions that determine optimal cutoff
combined <- function(model,data,n,reps,estimator,MAD,scale){

  #Use max sample size of 5000
  n <- min(n,5000)

  #apply function for simulated misspecified data
  misspec_fit <- miss_fit(model,data,n,reps,estimator,MAD,scale)

  #apply function for simulatied correct data
  true_fit <- true_fit(model,data,n,reps,estimator,MAD,scale)

  #Produce final table by level
  Table <- purrr::map(misspec_fit$misspec_fit_sum,~cbind(.,true_fit$true_fit_sum))

  #set up list of outcomes
  out<-list()

  #Table of fit indices
  out$Table<-Table
  #Table of discrepancies
  out$Q<-misspec_fit$all_Q
  #table of true generated data (for plotting distributions)
  out$D<-true_fit$all_D
  #table of individual simulated data (for plotting distributions)
  out$Data<-list(true_fit$all_D,misspec_fit$Data)

  return(out)
}

####################
#DDDFI function
####################
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

  if ( !(scale%in% c("normal", "nonnormal", "categorical"))){
    stop("dynamic Error: The scale of the observed variables must be 'normal','nonnormal', 'likert', or 'categorical'")}

  if (scale %in% "nonnormal"){
    warning("dynamic Warning:

            Computational times are longer for scale='nonnormal' due to missing data and procedures to generate non-normal data.

            This may take a few minutes.

            (e.g., 20 variables with N=1000 averages ~5-10 min, depending on the processor speed).", immediate.=T)
  }

  if (scale %in% "categorical"){
    warning("dynamic Warning:

            Computational times are longer for scale='categorical' due increased demand for categorical models.

            This may take a few minutes.

            (e.g., 20 variables with N=1000 averages ~5-10 min, depending on the processor speed).", immediate.=T)
  }

  #If model statement is manually entered,
  if(manual){

    # error for using manual=T with a lavaan object input
    tryCatch(cleanmodel(model),
             error=function(err5){
               if (grepl("no method for coercing this S4 class to a vector", err5)){
                 stop("dynamic Error: Did you accidentally include 'manual=TRUE' with a non-manually entered lavaan object?")
               }
             })

    n <- n
    model9 <- model

    #error for forgetting to include a sample size with manual=T
    tryCatch(defre(model9,n),
             error=function(err4){
               if (grepl("non-numeric matrix extent", err4)){
                 stop("dynamic Error: Did you forget to include a sample size with your manually entered model?")
               }
             })
  }

  #for model statment entered from a lavaan object,
  if(!manual){
    #error for entering a model statement manually but forgetting to include manual=T
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

  tryCatch(combined(model9,data,n,estimator, reps, MAD, scale, categorical),
           error=function(err3){
             if (grepl("is missing, with no default", err3)){
               stop("dynamic Error: Did you forget to include a dataset? 'scale' options other than 'normal' require a dataset to determine the number of categories/scale points.

                    If you do not have a dataset, scale='normal' will simulate complete data for all variables from normal distributions (cutoffs may be affected by missing data and deviations from normality) ")
             }
           })

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
                                                             ".96"="dotted"))+
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

  if(plot.dist & is.null(data)){
    stop("dynamic Error: Plots comparing generated data to the original data require a dataset to be provided.

           Please provide a original data with the 'data=' option or set the 'plot.dist=' option to FALSE
           ")
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
