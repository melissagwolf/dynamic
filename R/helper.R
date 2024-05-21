#' @importFrom dplyr filter select mutate full_join group_by add_tally
#' ungroup slice arrange summarise count as_tibble pull recode distinct
#' mutate_if slice_max slice_min
#' @importFrom tidyr pivot_longer unite nest extract
#' @importFrom lavaan lavaanify cfa fitMeasures
#' @importFrom simstandard sim_standardized
#' @importFrom purrr map map_dfr
#' @import ggplot2 GenOrd
#' @importFrom stats quantile qnorm cov2cor density rbeta rnorm runif
#' @importFrom semTools bsBootMiss
#' @importFrom MASS mvrnorm

############################################################################
########################### Multiple Functions #############################
############################################################################

value1<-density<-disc<-l<-value<-xx<-rand<-..density..<-NULL

#### Function to create model statement without numbers from user model (for input) ####
#Copy from OG
  cleanmodel <- function(model){

    suppressMessages(model %>%
                       lavaan::lavaanify(fixed.x = FALSE) %>%
                       dplyr::filter(lhs != rhs) %>%
                       dplyr::filter(op != "~1") %>%
                       dplyr::filter(op != "|") %>%
                       dplyr::group_by(lhs,op) %>%
                       dplyr::reframe(rhs = paste(rhs, collapse = " + ")) %>%
                       dplyr::arrange(dplyr::desc(op)) %>%
                       tidyr::unite("l", lhs, op, rhs, sep = " ") %>%
                       dplyr::pull(l))
  }

#### Function for Number of Factors ####
#Copy from OG

number_factor <- function(model){

  #prep the model
  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(lhs != rhs)

  #isolate factors
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #identify number of factors in model
  num_factors <- base::nrow(factors)

  return(num_factors)
}

#### Function for Isolating First-Order Factors ####

iso_first <- function(model){

  #isolate higher-order factors
  factHigh <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(op=="=~") %>%
    dplyr::filter(lhs %in% rhs)

  #if only up to second order, isolates first-order factors
  iso1 <- factHigh %>%
    dplyr::select(lhs) %>%
    base::unique() %>%
    base::unlist()

  #if more than two orders, isolates first-order factors
  iso2 <- factHigh %>%
    dplyr::filter(lhs %in% rhs) %>%
    dplyr::select(lhs) %>%
    base::unique() %>%
    base::unlist()

  #save object with first-order factors
  if(length(iso2) == 0) {factFirst <- iso1} else {factFirst <- iso2}

  return(factFirst)
}

#### Function for Number of First-Order Factors ####

number_factor_first <- function(model){

  #identify number of factors in model
  num_factors <- base::length(iso_first(model))

  return(num_factors)
}

#Did they enter unstandardized loadings?  Aka, do they have any loadings = 1?
#Copy from OG
#Used for error message

unstandardized <- function(model){

  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(lhs != rhs)

  one_plus <- lav_file %>%
    dplyr::filter(ustart >= 1 | ustart <= -1) %>%
    base::nrow()

  return(one_plus)
}

#### Function to calculate degrees of freedom ####
# Used for error message

defre <- function(model,n){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Rename
  true_dgm <- model

  #Run one simulation
  dat <- simstandard::sim_standardized(true_dgm,n=n,latent=FALSE,errors=FALSE)
  fit <- lavaan::cfa(model=mod,data=dat,std.lv=TRUE)

  #Number of freely estimated paths
  paths <- base::max(lavaan::parTable(fit)$free)

  #Number of unique values in input matrix
  parms <- base::nrow(lavaan::lavInspect(fit,"std.lv")$theta)
  tot.parms <- (parms*(1+parms))/2

  #Subtract
  return(tot.parms-paths)
}

#SPECIFIC TO R PACKAGE

cfa_n <- function(model){

  #Extract n from lavaan object
  #Warning message to hide warning when someone enters a non-lavaan object (error message will display instead)
  #Only need it here because this is the first argument in cfaHB and cfaOne
  n <- base::unlist(model@SampleStats@nobs)
  return(n)
}

##### Extract model statement from lavaan object #####
#SPECIFIC TO R PACKAGE

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
                               dplyr::reframe(rhs=paste(est.std,"*",rhs,collapse=" + ")) %>%
                               dplyr::arrange(desc(op)) %>%
                               tidyr::unite("mod",lhs,op,rhs,sep="") %>%
                               dplyr::pull(mod))

  #Collapse into one string because my other functions expect that
  mod <- base::paste(ss_mod, sep="", collapse="\n")

  return(mod)

}

#### Function to clean dataset of fit indices (now exportable as of 1.26.22)

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


############################################################################
################################# cfaOne ###################################
############################################################################



### One-factor: Function to see which items are available ###

one_num <- function(model){

  #Rename (just to be consistent with shiny app)
  Mod_C <- model

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(Mod_C, fixed.x=FALSE) %>%
    dplyr::filter(lhs != rhs)

  #identify all factor names
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #Identify any items that already have an error covariance
  items_covariance <- factors %>%
    dplyr::mutate(type="Factor") %>%
    dplyr::full_join(lav_file, by = "lhs") %>%
    dplyr::select(-type,type) %>%
    dplyr::select(lhs,op,rhs,type) %>%
    dplyr::filter(op=="=~" | is.na(type)) %>%
    dplyr::filter(is.na(type)) %>%
    dplyr::select(-type) %>%
    tidyr::pivot_longer(-op,names_to = "test", values_to = "rhs") %>%
    dplyr::select(-op,-test) %>%
    dplyr::mutate(lhs=NA,op=NA,ustart=NA)

  #Isolate the items that do not already have an error covariance
  solo_items <- lav_file %>%
    dplyr::select(lhs,op,rhs,ustart) %>%
    base::rbind(items_covariance) %>%
    dplyr::filter(op=="=~"|is.na(op)) %>%
    dplyr::group_by(rhs) %>%
    dplyr::add_tally() %>%
    dplyr::filter(n==1) %>%
    dplyr::ungroup() %>%
    arrange(abs(ustart))

  return(solo_items)
}

#### One-Factor: Function to create misspecification statement ####

one_add <- function(model){

  #Read in available items
  itemoptions <- one_num(model)

  #Count number of available items
  num_i <- base::nrow(itemoptions)

  #Select items for misspecification depending on number of available items
  if(num_i==4){
    num_m <- itemoptions %>%
      dplyr::slice(1:2)
  }else if(num_i==5){
    num_m <- itemoptions %>%
      dplyr::slice(1:4)
  }else{
    num_m <- itemoptions %>%
      dplyr::slice(1:(floor(num_i/2)*2))
  }

  #Identifiers to separate odds and even rows
  evenindex <- base::seq(2,base::nrow(num_m),2)
  oddindex <- base::seq(1,base::nrow(num_m),2)

  #Separate
  left <- num_m[evenindex,]
  right <- num_m[oddindex,] %>%
    `colnames<-`(c("lhs_1","op_1","rhs_1","ustart_1","n_1"))

  #Create misspecification statements
  Residual_Correlation <- base::cbind(left,right) %>%
    dplyr::mutate(cor=.3,
                  opp="~~",
                  star="*") %>%
    tidyr::unite(V1,c("rhs","opp","cor","star","rhs_1"),sep=" ") %>%
    dplyr::select(V1)

  return(Residual_Correlation)
}

#### One-factor: Function to create Misspecified DGM ####

DGM_one <- function(model){

  #Count number of available items for number of misspecifications
  num_m<- base::nrow(one_num(model))

  #Figure out number of levels given number of available items
  if(num_m==4){
    L1 <- 1
    levels <- L1
  }else if(num_m==5){
    L1 <- 1
    L2 <- 2
    levels <- base::rbind(L1,L2)
  }else{
    L3 <- base::floor(num_m/2)
    L2 <- base::floor((2*L3)/3)
    L1 <- base::floor(L3/3)
    levels <- base::rbind(L1,L2,L3)
  }

  #Read in misspecifications
  mod <- one_add(model)

  #Get parameters for true dgm
  Mods <- model
  #Mod_C <- base::as.character(Mods$V1)

  #single_mod <- base::lapply(levels, function(x) base::rbind(Mod_C,mod[base::seq(x), ,drop = FALSE]) %>%
  #                       base::data.frame() %>%
  #                       dplyr::pull(V1))
  #This made you miserable. Shiny was struggling with \n at the end of strings here, for some reason.

  #Create a list for every row in the mod object (misspecifications)
  #For each element, bind the misspecification to the OG model statement sequentially
  #Turn it into a dataframe and extract
  single_mod <- base::lapply(levels, function(x) base::rbind(Mods,mod[base::seq(x), ,drop = FALSE]) %>%
                               base::data.frame() %>%
                               dplyr::pull(V1))

  return(single_mod)

}

### One-factor: Simulate fit indices for misspecified model for all levels ###

one_fit <- function(model,n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm
  misspec_dgm <- DGM_one(model)

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set seed
  set.seed(649364)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Simulate one large dataset for each misspecification (use map to apply across each
  #element (set of misspecifications) in the list)
  all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*r,
                                                                            latent=FALSE,errors=FALSE))

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_misspec <- purrr::map(all_data_misspec,~base::cbind(.,rep_id_misspec))

  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~dplyr::group_by(.,rep_id_misspec) %>%
                               tidyr::nest())

  #Grab data level of the list
  data <- purrr::map(misspec_data,2)

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa for each element in the list
  misspec_cfa <- purrr::map(data, function(x) purrr::map(x, function(y) lavaan::cfa(model=mod,
                                                                                    estimator=estimator,
                                                                                    data=y,
                                                                                    std.lv=TRUE,
                                                                                    se=se,
                                                                                    check.gradient=FALSE,
                                                                                    check.post=FALSE,
                                                                                    check.vcov=FALSE,
                                                                                    control=list(rel.tol=.001))))

  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, ind)) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))

  set.seed(NULL)

  return(misspec_fit_sum)

}

#### One_Factor: Function to create True DGM (aka, just the model the user read in) ####

true_fit_one <- function(model,n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  true_dgm <- model

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set Seed
  set.seed(326267)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Simulate one large dataset
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n*r,
                                                 latent = FALSE,
                                                 errors = FALSE)

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_true <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)

  #Group and list
  true_data <- dat_rep_true %>%
    dplyr::group_by(rep_id_true) %>%
    tidyr::nest() %>%
    base::as.list()

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa
  true_cfa <- purrr::map(true_data$data,function(x) lavaan::cfa(model=mod,
                                                                estimator=estimator,
                                                                data=x,
                                                                std.lv=TRUE,
                                                                se=se,
                                                                check.gradient=FALSE,
                                                                check.post=FALSE,
                                                                check.vcov=FALSE,
                                                                control=list(rel.tol=.001)))

  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., ind)) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#### One-Factor: Function to combine both model fit stats for all levels into one dataframe ####

one_df <- function(model,n,estimator,reps){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- one_fit(model,n,estimator,reps)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_one(model,n,estimator,reps)

  #Produce final table by level
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}


############################################################################
################################# cfaHB ####################################
############################################################################


### Multi-factor: Function to see which items are available ###

multi_num_HB <- function(model){

  #Rename (just to be consistent with shiny app)
  Mod_C <- model

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(Mod_C, fixed.x=FALSE) %>%
    dplyr::filter(lhs != rhs)

  #identify all factor names
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #Identify number of items per factor
  num_items <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::group_by(lhs) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    base::as.data.frame() %>%
    `colnames<-`(c("lhs","Original"))

  #Identify any items that already have an error covariance
  items_covariance <- factors %>%
    dplyr::mutate(type="Factor") %>%
    dplyr::full_join(lav_file, by = "lhs") %>%
    dplyr::select(-type,type) %>%
    dplyr::select(lhs,op,rhs,type) %>%
    dplyr::filter(op=="=~" | is.na(type)) %>%
    dplyr::filter(is.na(type)) %>%
    dplyr::select(-type) %>%
    tidyr::pivot_longer(-op,names_to = "test", values_to = "rhs") %>%
    dplyr::select(-op,-test) %>%
    dplyr::mutate(lhs=NA,op=NA,ustart=NA)

  #Isolate the items that do not already have an error covariance or cross-loading
  solo_items <- lav_file %>%
    dplyr::select(lhs,op,rhs,ustart) %>%
    base::rbind(items_covariance) %>%
    dplyr::filter(op=="=~"|is.na(op)) %>%
    dplyr::group_by(rhs) %>%
    dplyr::add_tally() %>%
    dplyr::filter(n==1) %>%
    dplyr::ungroup()

  #Count number of items remaining per factor
  remaining <- solo_items %>%
    dplyr::group_by(lhs) %>%
    dplyr::select(-n) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    dplyr::full_join(num_items,by="lhs") %>%
    base::as.data.frame() %>%
    `colnames<-`(c("lhs","Remaining","Original"))

  #Add in factor loadings, group by number of items per factor (>2 or 2)
  #And sort factor loadings magnitude within group
  itemoptions <- solo_items %>%
    dplyr::full_join(remaining,by="lhs") %>%
    dplyr::mutate(priority=ifelse(Original>2 & Remaining !="NA","Three","Two")) %>%
    dplyr::group_by(priority) %>%
    dplyr::arrange(abs(ustart), .by_group=TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::select(lhs,rhs,ustart,priority) %>%
    base::as.data.frame() %>%
    dplyr::as_tibble() %>%
    `colnames<-`(c("lhs","Item","Loading","Priority"))

  return(itemoptions)
}

#### Multi-Factor: Function to identify available items and loading magnitude ####

multi_add_HB <- function(model){

  #read in the model
  Mod_C <- model

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(Mod_C, fixed.x=FALSE) %>%
    dplyr::filter(lhs != rhs)

  #read in number of factors
  num_fact <- number_factor(model)

  #read in viable items from each factor
  itemoptions <- multi_num_HB(model)

  #select lowest loading from each factor, in order of magnitude
  crosses <- itemoptions %>%
    dplyr::group_by(lhs) %>%
    dplyr::slice_min(base::abs(Loading)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Priority,base::abs(Loading)) %>%
    dplyr::slice(1:(num_fact-1))

  #identify all factor names (again)
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #Compute Coefficient H for each factor
  suppressMessages(Coef_H <- lavaan::lavaanify(Mod_C, fixed.x = FALSE) %>%
                     dplyr::filter(lhs != rhs) %>%
                     dplyr::filter(op == "=~") %>%
                     dplyr::mutate(L_Sq=ustart^2) %>%
                     dplyr::mutate(E_Var=1-L_Sq) %>%
                     dplyr::mutate(Div=L_Sq/E_Var) %>%
                     dplyr::group_by(lhs) %>%
                     dplyr::reframe(Sum=sum(Div)) %>%
                     dplyr::mutate(rand=runif(n=nrow(factors),0,.001)) %>%
                     dplyr::mutate(H=((1+(Sum^-1))^-1)+rand) %>%
                     dplyr::select(-Sum,-rand) %>%
                     dplyr::arrange(-H) %>%
                     `colnames<-`(c("rhs","H")))

  #isolate factors and factor correlations
  factcor1 <- factors %>%
    dplyr::mutate(type="Factor") %>%
    dplyr::full_join(lav_file, by = "lhs") %>%
    dplyr::mutate(type=dplyr::recode(type, .missing ="Error Correlation")) %>%
    dplyr::select(lhs,op,rhs,ustart,type) %>%
    dplyr::filter(op=="~~" & type=="Factor")

  #flip in reverse so we get a list of all factors in one column
  factcor2 <- factors %>%
    dplyr::mutate(type="Factor") %>%
    dplyr::full_join(lav_file, by = "lhs") %>%
    dplyr::select(lhs,op,rhs,ustart,type) %>%
    dplyr::filter(op=="~~" & type=="Factor") %>%
    `colnames<-`(c("rhs","op","lhs","ustart","type")) %>%
    dplyr::select(lhs,op,rhs,ustart,type)

  #Isolate items
  dup1 <- factcor1 %>%
    dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>%
    dplyr::full_join(crosses,by="lhs") %>%
    dplyr::full_join(Coef_H,by="rhs") %>%
    dplyr::filter(Item != "NA") %>%
    dplyr::arrange(abs(Loading))

  #Run twice for cleaning later
  dup2 <- dup1

  #Manipulate to create model statement
  #Need to add factor correlation statement where lowest comes first
  #So that we can remove from contention once a factor correlation is added
  setup <- base::rbind(dup1,dup2) %>%
    dplyr::mutate(lhs_1=lhs,
                  rhs_1=rhs,
                  f_min=base::pmin(lhs_1,rhs_1),
                  f_max=base::pmax(lhs_1,rhs_1)) %>%
    tidyr::unite(facts,c("f_min","f_max")) %>%
    dplyr::select(-lhs_1,-rhs_1) %>%
    dplyr::distinct(lhs,op,rhs,ustart,type,Item,Loading,Priority,H,.keep_all = TRUE)


  #Rename for iteration
  setup_copy <- setup

  #Create empty dataframe
  cleaned <- base::data.frame(base::matrix(nrow=0,ncol=10)) %>%
    `colnames<-`(names(setup)) %>%
    dplyr::mutate_if(is.logical, as.character)

  #Cycle through to grab F-1 misspecifications
  #Select the highest H for first item I (for crossloading)
  #Use anti_join to remove that factor correlation from the list for the next item
  for (i in unique(setup_copy$Item)){
    cleaned[i,] <- setup_copy %>%
      dplyr::filter(Item==i) %>%
      dplyr::slice_max(H)
    setup_copy <- dplyr::anti_join(setup_copy,cleaned,by="facts")
  }

  #Prep dtaframe for model statement
  modinfo <- cleaned %>%
    dplyr::mutate(operator="=~",
                  H=base::as.numeric(H),
                  Loading=base::as.numeric(Loading),
                  ustart=base::as.numeric(ustart)) %>%
    dplyr::arrange(Priority,Loading,-H)

  #Compute maximum allowable cross loading value
  Cross_Loading <- modinfo %>%
    dplyr::mutate(F1=ustart,
                  F1_Sq=F1^2,
                  L1=Loading,
                  L1_Sq=L1^2,
                  E=1-L1_Sq) %>%
    dplyr::mutate(MaxAllow=((base::sqrt(((L1_Sq*F1_Sq)+E))-(abs(L1*F1)))*.95),
                  MaxAllow2=base::round(MaxAllow,digits=4),
                  Final_Loading=base::pmin(abs(Loading),abs(MaxAllow2)),
                  times="*") %>%
    dplyr::select(rhs,operator,Final_Loading,times,Item) %>%
    tidyr::unite("V1",sep=" ")

  #return value to append to model statement
  return(Cross_Loading)
}


#### Multi-factor: Function to create Misspecified DGM given the number of factors ####

DGM_Multi_HB <- function(model){

  mod <- multi_add_HB(model)

  #Get parameters for true dgm
  Mods <- model
  #Mod_C <- base::as.character(Mods$V1)

  #multi_mod <- lapply(mod, function(x) rbind(Mod_C,mod[seq(x), ,drop = FALSE]) %>%
  #                      data.frame() %>%
  #                      pull(V1))

  #This made you miserable. Shiny/R was struggling with \n at the end of strings here, for some reason.

  #Create a list for every row in the mod object (misspecifications)
  #For each element, bind the misspecification to the OG model statement sequentially
  #Turn it into a dataframe and extract
  multi_mod <- lapply(seq(nrow(mod)), function(x) rbind(Mods,mod[seq(x), ,drop = FALSE]) %>%
                        base::data.frame() %>%
                        dplyr::pull(V1))

  return(multi_mod)

}

### Multi-factor: Simulate fit indices for misspecified model for all levels ###

multi_fit_HB <- function(model,n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm (this is a list)
  misspec_dgm <- DGM_Multi_HB(model)

  #Use max sample size of 2000
  n <- min(n,2000)

  #Set seed
  set.seed(269854)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Simulate one large dataset for each misspecification (use map to apply across each
  #element (set of misspecifications) in the list)
  all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*r,
                                                                            latent=FALSE,errors=FALSE))

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- rep(1:r,n)

  #Combine indicator with dataset for each element in list
  dat_rep_misspec <- purrr::map(all_data_misspec,~cbind(.,rep_id_misspec))

  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~group_by(.,rep_id_misspec) %>%
                               tidyr::nest())

  #Grab data level of the list
  data <- purrr::map(misspec_data,2)

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa for each element in the list
  misspec_cfa <- purrr::map(data, function(x) purrr::map(x, function(y) lavaan::cfa(model=mod,
                                                                                    estimator=estimator,
                                                                                    data=y,
                                                                                    std.lv=TRUE,
                                                                                    se=se,
                                                                                    check.gradient=FALSE,
                                                                                    check.post=FALSE,
                                                                                    check.vcov=FALSE,
                                                                                    control=list(rel.tol=.001))))

  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, ind)) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))
  set.seed(NULL)

  return(misspec_fit_sum)

}

#### Multi_Factor: Function to create True DGM (aka, just the model the user read in) ####

true_fit_HB <- function(model,n,estimator,reps){

  #Can make this faster by only doing it once
  #Would need to change table. Not sure what would happen to plot.
  #Already did this

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for true DGM
  true_dgm <- model

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Set Seed
  set.seed(267326)

  #Simulate one large dataset
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n*r,
                                                 latent = FALSE,
                                                 errors = FALSE)

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_true <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)

  #Group and list
  true_data <- dat_rep_true %>%
    dplyr::group_by(rep_id_true) %>%
    tidyr::nest() %>%
    base::as.list()

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa
  true_cfa <- purrr::map(true_data$data,function(x) lavaan::cfa(model=mod,
                                                                estimator=estimator,
                                                                data=x,
                                                                std.lv=TRUE,
                                                                se=se,
                                                                check.gradient=FALSE,
                                                                check.post=FALSE,
                                                                check.vcov=FALSE,
                                                                control=list(rel.tol=.001)))

  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., ind)) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#### Multi-Factor: Function to combine both model fit stats for all levels into one dataframe ####

multi_df_HB <- function(model,n,estimator,reps){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- multi_fit_HB(model,n,estimator,reps)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_HB(model,n,estimator,reps)

  #Produce final table of fit indices for each level (as a list)
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}

############################################################################
################################ equivTest #################################
############################################################################

##### NCP Chi-Sq #####

equiv_ncp_chi2 <- function(alpha,T_ml,df){

  z=stats::qnorm(1-alpha)
  z2=z*z
  z3=z2*z
  z4=z3*z
  z5=z4*z
  sig2=2*(2*T_ml-df+2)
  sig=sqrt(sig2)
  sig3=sig*sig2
  sig4=sig2*sig2
  sig5=sig4*sig
  sig6=sig2*sig4

  delta=T_ml-df+2+sig*
    (
      z+(z2-1)/sig-z/sig2 + 2*(df-1)*(z2-1)/(3*sig3)
      +( -(df-1)*(4*z3-z)/6+(df-2)*z/2 )/sig4
      +4*(df-1)*(3*z4+2*z2-11)/(15*sig5)
      +(
        -(df-1)*(96*z5+164*z3-767*z)/90-4*(df-1)*(df-2)*(2*z3-5*z)/9
        +(df-2)*z/2
      )/sig6
    )

  delta=max(delta,0)

  return(delta)
}

######## Compute Values #########

equiv_cutoffs <- function(p,T_ml,df,T_mli,n){

  #Set parms
  df_i <- p*(p-1)/2
  alpha <- .05

  #T-size RMSEA#;
  delta_t_r <- equiv_ncp_chi2(alpha,T_ml,df)
  RMSEA_t <- sqrt(delta_t_r/(df*n))

  #T-size CFI

  delta_t_c <- equiv_ncp_chi2(alpha/2, T_ml,df)
  delta_it <- equiv_ncp_chi2(1-alpha/2, T_mli,df_i)
  CFI_t <- 1-max(delta_t_c,0)/max(delta_t_c,delta_it,0)

  #Recalculate Bins based on Model Characteristics - RMSEA#

  RMSEA_e01=exp(
    1.34863-.51999*log(df)+.01925*log(df)*log(df)-.59811*log(n)+.00902*sqrt(n)+.01796*log(df)*log(n))


  RMSEA_e05=exp(2.06034-.62974*log(df)+.02512*log(df)*log(df)-.98388*log(n)
                +.05442*log(n)*log(n)-.00005188*n+.05260*log(df)*log(n))


  RMSEA_e08=exp(2.84129-.54809*log(df)+.02296*log(df)*log(df)-.76005*log(n)
                +.10229*log(n)*log(n)-1.11167*(n^.2)+.04845*log(df)*log(n))


  RMSEA_e10=exp(2.36352-.49440*log(df)+.02131*log(df)*log(df)-.64445*log(n)
                +.09043*log(n)*log(n)-1.01634*(n^.2)+.04422*log(df)*log(n))

  ## Recalculate - CFI

  CFI_e99=1-exp(
    4.67603-.50827*log(df)+.87087*(df^(1/5))-.59613*((df_i)^(1/5))-1.89602*log(n)
    + .10190*((log(n))^2)+ .03729*log(df)*log(n)
  );

  #corresponding to R-square=.9836;

  CFI_e95=1-exp(
    4.12132-.46285*log(df)+.52478*(df^(1/5))-.31832*((df_i)^(1/5))-1.74422*log(n)
    +.13042*((log(n))^2)-.02360*(n^(1/2))+.04215*log(df)*log(n)
  );

  #corresponding to R-square=.9748;

  CFI_e92=1-exp(
    6.31234-.41762*log(df)+.01554*((log(df))^2)-.00563*((log(df_i))^2)-1.30229*log(n)
    +.19999*((log(n))^2)-2.17429*(n^(1/5))+.05342*log(df)*log(n)-.01520*log(df_i)*log(n)
  );

  #corresponding to R-square=.9724

  CFI_e90=1-exp(
    5.96633-.40425*log(df)+.01384*((log(df))^2)-.00411*((log(df_i))^2)-1.20242*log(n)
    +.18763*((log(n))^2)-2.06704*(n^(1/5))+.05245*log(df)*log(n)-.01533*log(df_i)*log(n)
  );

  ## Create bins

  cutoff_rmsea <- cbind(RMSEA_e01, RMSEA_e05, RMSEA_e08, RMSEA_e10, RMSEA_t)
  cutoff_cfi <- cbind(CFI_e90, CFI_e92, CFI_e95, CFI_e99, CFI_t)
  cutoff_combo <- rbind(cutoff_rmsea,cutoff_cfi)
  cutoff_3 <- round(cutoff_combo,3)
  colnames(cutoff_3) <- c("Cut_1","Cut_2","Cut_3","Cut_4","T")

  return(cutoff_3)
}

####### Lavaan extraction #######

equiv_n <- function(obj){
  n <- base::unlist(obj@SampleStats@nobs)
  return(n)
}

equiv_p <- function(obj){
  p <- base::unlist(obj@Model@nvar)
  return(p)
}

equiv_T_ml <- function(obj){
  #Warning message to hide warning when someone enters a non-lavaan object (error message will display instead)
  #Only need it here because this is the first argument in equivTest
  T_ml <- base::unlist(obj@test[["standard"]][["stat"]])
  return(T_ml)
}

equiv_T_mli <- function(obj){
  T_mli <- base::unlist(obj@baseline[["test"]][["standard"]][["stat"]])
  return(T_mli)
}

equiv_df <- function(obj){
  df <- base::unlist(obj@test[["standard"]][["df"]])
  return(df)
}

############################################################################
################################ exactFit ##################################
############################################################################

exact_fit_dat <- function(model,n,reps){
  mod <- cleanmodel(model)

  true_dgm <- model

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set Seed
  set.seed(267326)

  #Number of reps
  r <- reps

  #Simulate one large dataset
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n*r,
                                                 latent = FALSE,
                                                 errors = FALSE)

  #Create indicator to split into r datasets for r reps
  rep_id_true <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)

  #Group and list
  true_data <- dat_rep_true %>%
    dplyr::group_by(rep_id_true) %>%
    tidyr::nest() %>%
    base::as.list()

  #Run 500 cfa
  true_sem <- purrr::map(true_data$data,~base::withCallingHandlers(lavaan::sem(model = mod, data=., std.lv=TRUE),
                                                                   check.gradient=FALSE,
                                                                   check.post=FALSE,
                                                                   check.vcov=FALSE))

  #Extract fit stats from each rep (list) into a data frame and clean
  true_fit_sum <- purrr::map_dfr(true_sem,~lavaan::fitMeasures(., c("chisq", "df","srmr","rmsea","cfi")))

  set.seed(NULL)

  return(true_fit_sum)
}

############################################################################
################################# hier1HB ##################################
############################################################################

#### Function to Create lav_file with Only First-Order Factors ####

lav_file_first <- function(model){

  #isolate names of first-order factors
  factFirst <- iso_first(model)

  #isolate names of second-order factors
  factSecond <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(op=="=~") %>%
    dplyr::filter(rhs %in% factFirst) %>%
    dplyr::select(lhs) %>%
    base::unique() %>%
    base::unlist()

  #separate out factor correlation rows
  corr_lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(op=="~~") %>%
    dplyr::filter(rhs != lhs) %>%
    dplyr::filter(lhs %in% factFirst) %>%
    dplyr::filter(rhs %in% factFirst)

  #isolate first-order factor rows
  fact_lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(op=="=~") %>%
    dplyr::filter(lhs %in% factFirst)

  #isolate second-order factor rows
  rowsSecond <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(op=="=~") %>%
    dplyr::filter(lhs %in% factSecond)

  #isolate rows with items that cross-load onto first- and second-order factors
  cross_lav_file <- rowsSecond %>% filter(rowsSecond$rhs %in% fact_lav_file$rhs)

  #prep the model
  lav_file <- dplyr::bind_rows(cross_lav_file, fact_lav_file, corr_lav_file)

  return(lav_file)
}

#### Function for Model-Implied First-Order Factor Correlations  ####

MI_corr_first <- function(model){

  #isolate first-order factors
  factFirst <- iso_first(model)

  #all possible pairs of second-order factors
  pairs <- base::t(utils::combn(factFirst, 2))

  #initialize data frame for model-implied correlations
  df <- dplyr::tibble(lhs = pairs[,1], rhs = pairs[,2]) %>%
    dplyr::mutate(op = "~~", type = "Factor") %>%
    dplyr::relocate(op, .after = lhs)

  #compute all model-implied correlations
  impliedCorr <- simstandard::sim_standardized_matrices(model)
  mat <- impliedCorr$Correlations$R

  #collect model-implied correlations for second-order factors
  factcor <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(ustart = mat[lhs,rhs]) %>%
    dplyr::relocate(ustart, .before = type)

  #create df to mimic factcor1 in multi_add_HB
  factcor1 <- factcor %>%
    dplyr::relocate(ustart, .before = type)

  return(factcor1)
}

### Hierarchical: Function to see which items are available ###

multi_num_hier <- function(model){

  #Rename (just to be consistent with shiny app)
  Mod_C <- model

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lav_file_first(Mod_C)

  #identify all factor names
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #Identify number of items per factor
  num_items <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::group_by(lhs) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    base::as.data.frame() %>%
    `colnames<-`(c("lhs","Original"))

  #Identify any items that already have an error covariance
  items_covariance <- factors %>%
    dplyr::mutate(type="Factor") %>%
    dplyr::full_join(lav_file, by = "lhs") %>%
    dplyr::select(-type,type) %>%
    dplyr::select(lhs,op,rhs,type) %>%
    dplyr::filter(op=="=~" | is.na(type)) %>%
    dplyr::filter(is.na(type)) %>%
    dplyr::select(-type) %>%
    tidyr::pivot_longer(-op,names_to = "test", values_to = "rhs") %>%
    dplyr::select(-op,-test) %>%
    dplyr::mutate(lhs=NA,op=NA,ustart=NA)

  #Isolate the items that do not already have an error covariance or cross-loading
  solo_items <- lav_file %>%
    dplyr::select(lhs,op,rhs,ustart) %>%
    base::rbind(items_covariance) %>%
    dplyr::filter(op=="=~"|is.na(op)) %>%
    dplyr::group_by(rhs) %>%
    dplyr::add_tally() %>%
    dplyr::filter(n==1) %>%
    dplyr::ungroup()

  #Count number of items remaining per factor
  remaining <- solo_items %>%
    dplyr::group_by(lhs) %>%
    dplyr::select(-n) %>%
    dplyr::count() %>%
    dplyr::ungroup() %>%
    dplyr::full_join(num_items,by="lhs") %>%
    base::as.data.frame() %>%
    `colnames<-`(c("lhs","Remaining","Original"))

  #Add in factor loadings, group by number of items per factor (>2 or 2)
  #And sort factor loadings magnitude within group
  itemoptions <- solo_items %>%
    dplyr::full_join(remaining,by="lhs") %>%
    dplyr::mutate(priority=ifelse(Original>2 & Remaining !="NA","Three","Two")) %>%
    dplyr::group_by(priority) %>%
    dplyr::arrange(abs(ustart), .by_group=TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::select(lhs,rhs,ustart,priority) %>%
    base::as.data.frame() %>%
    dplyr::as_tibble() %>%
    `colnames<-`(c("lhs","Item","Loading","Priority"))

  return(itemoptions)
}

#### Hierarchical: Function to identify available items and loading magnitude ####

multi_add_hier <- function(model){

  #read in the model
  Mod_C <- model

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lav_file_first(Mod_C)

  #read in number of factors
  num_fact <- number_factor_first(model)

  #read in viable items from each factor
  itemoptions <- multi_num_hier(model)

  #select lowest loading from each factor, in order of magnitude
  crosses <- itemoptions %>%
    dplyr::group_by(lhs) %>%
    dplyr::slice_min(base::abs(Loading)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Priority,base::abs(Loading)) %>%
    dplyr::slice(1:(num_fact-1))

  #identify all factor names (again)
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #Compute Coefficient H for each factor
  suppressMessages(Coef_H <- lav_file %>%
                     dplyr::filter(op == "=~") %>%
                     dplyr::mutate(L_Sq=ustart^2) %>%
                     dplyr::mutate(E_Var=1-L_Sq) %>%
                     dplyr::mutate(Div=L_Sq/E_Var) %>%
                     dplyr::group_by(lhs) %>%
                     dplyr::summarise(Sum=sum(Div)) %>%
                     dplyr::mutate(H=((1+(Sum^-1))^-1)) %>%
                     dplyr::select(-Sum) %>%
                     dplyr::arrange(-H) %>%
                     `colnames<-`(c("rhs","H")))

  #isolate factors and factor correlations
  factcor1 <- MI_corr_first(model)

  #flip in reverse so we get a list of all factors in one column
  factcor2 <- factcor1
  factcor2 <- dplyr::rename(factcor2, lhs = rhs, rhs = lhs) %>%
    dplyr::relocate(rhs, .after = op) %>%
    dplyr::relocate(lhs, .before = op)

  #Isolate items
  dup1 <- factcor1 %>%
    dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>%
    dplyr::full_join(crosses,by="lhs") %>%
    dplyr::full_join(Coef_H,by="rhs") %>%
    dplyr::filter(Item != "NA") %>%
    dplyr::arrange(abs(Loading))

  #Run twice for cleaning later
  dup2 <- dup1

  #Manipulate to create model statement
  #Need to add factor correlation statement where lowest comes first
  #So that we can remove from contention once a factor correlation is added
  setup <- base::rbind(dup1,dup2) %>%
    dplyr::mutate(lhs_1=lhs,
                  rhs_1=rhs) %>%
    dplyr::mutate(f_min=base::pmin(lhs_1,rhs_1),
                  f_max=base::pmax(lhs_1,rhs_1)) %>%
    tidyr::unite(facts,c("f_min","f_max")) %>%
    dplyr::select(-lhs_1,-rhs_1) %>%
    dplyr::distinct(lhs,op,rhs,ustart,type,Item,Loading,Priority,H,.keep_all = TRUE)

  #Rename for iteration
  setup_copy <- setup

  #Create empty dataframe
  cleaned <- base::data.frame(base::matrix(nrow=0,ncol=10)) %>%
    `colnames<-`(names(setup)) %>%
    dplyr::mutate_if(is.logical, as.character)

  #Cycle through to grab F-1 misspecifications
  #Select the highest H for first item I (for crossloading)
  #Use anti_join to remove that factor correlation from the list for the next item
  for (i in unique(setup_copy$Item)){
    cleaned[i,] <- setup_copy %>%
      dplyr::filter(Item==i) %>%
      dplyr::slice_max(H)
    setup_copy <- dplyr::anti_join(setup_copy,cleaned,by="facts")
  }

  #Prep dtaframe for model statement
  modinfo <- cleaned %>%
    dplyr::mutate(operator="=~",
                  H=base::as.numeric(H),
                  Loading=base::as.numeric(Loading),
                  ustart=base::as.numeric(ustart)) %>%
    dplyr::arrange(Priority,Loading,-H)

  #Compute maximum allowable cross loading value
  Cross_Loading <- modinfo %>%
    dplyr::mutate(F1=ustart,
                  F1_Sq=F1^2,
                  L1=Loading,
                  L1_Sq=L1^2,
                  E=1-L1_Sq) %>%
    dplyr::mutate(MaxAllow=((base::sqrt(((L1_Sq*F1_Sq)+E))-(abs(L1*F1)))*.95),
                  MaxAllow2=base::round(MaxAllow,digits=4),
                  Final_Loading=base::pmin(abs(Loading),abs(MaxAllow2)),
                  times="*") %>%
    dplyr::select(rhs,operator,Final_Loading,times,Item) %>%
    tidyr::unite("V1",sep=" ")

  #return value to append to model statement
  return(Cross_Loading)
}

#### Hierarchical: Function to create Misspecified DGM given the number of factors ####

DGM_Multi_hier <- function(model){

  mod <- multi_add_hier(model)

  #Get parameters for true dgm
  Mods <- model
  #Mod_C <- base::as.character(Mods$V1)

  #multi_mod <- lapply(mod, function(x) rbind(Mod_C,mod[seq(x), ,drop = FALSE]) %>%
  #                      data.frame() %>%
  #                      pull(V1))

  #This made you miserable. Shiny/R was struggling with \n at the end of strings here, for some reason.

  #Create a list for every row in the mod object (misspecifications)
  #For each element, bind the misspecification to the OG model statement sequentially
  #Turn it into a dataframe and extract
  multi_mod <- lapply(seq(nrow(mod)), function(x) rbind(Mods,mod[seq(x), ,drop = FALSE]) %>%
                        base::data.frame() %>%
                        dplyr::pull(V1))

  return(multi_mod)

}

### Hierarchical: Simulate fit indices for misspecified model for all levels ###

multi_fit_hier <- function(model,n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm (this is a list)
  misspec_dgm <- DGM_Multi_hier(model)

  #Use max sample size of 2000
  n <- min(n,2000)

  #Set seed
  set.seed(269854)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Simulate one large dataset for each misspecification (use map to apply across each
  #element (set of misspecifications) in the list)
  all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*r,
                                                                            latent=FALSE,errors=FALSE))

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- rep(1:r,n)

  #Combine indicator with dataset for each element in list
  dat_rep_misspec <- purrr::map(all_data_misspec,~cbind(.,rep_id_misspec))

  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~group_by(.,rep_id_misspec) %>%
                               tidyr::nest())

  #Grab data level of the list
  data <- purrr::map(misspec_data,2)

  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa for each element in the list
  misspec_cfa <- purrr::map(data, function(x) purrr::map(x, function(y) lavaan::cfa(model=mod,
                                                                                    estimator=estimator,
                                                                                    data=y,
                                                                                    std.lv=TRUE,
                                                                                    se=se,
                                                                                    check.gradient=FALSE,
                                                                                    check.post=FALSE,
                                                                                    check.vcov=FALSE,
                                                                                    control=list(rel.tol=.001))))

  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, ind)) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))


  set.seed(NULL)

  return(misspec_fit_sum)

}

#### Hierarchical: Function to create True DGM (aka, just the model the user read in) ####

true_fit_hier <- function(model,n,estimator,reps){

  #Can make this faster by only doing it once
  #Would need to change table. Not sure what would happen to plot.
  #Already did this

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for true DGM
  true_dgm <- model

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Set Seed
  set.seed(267326)

  #Simulate one large dataset
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n*r,
                                                 latent = FALSE,
                                                 errors = FALSE)

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_true <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)

  #Group and list
  true_data <- dat_rep_true %>%
    dplyr::group_by(rep_id_true) %>%
    tidyr::nest() %>%
    base::as.list()


  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa
  true_cfa <- purrr::map(true_data$data,function(x) lavaan::cfa(model=mod,
                                                                estimator=estimator,
                                                                data=x,
                                                                std.lv=TRUE,
                                                                se=se,
                                                                check.gradient=FALSE,
                                                                check.post=FALSE,
                                                                check.vcov=FALSE,
                                                                control=list(rel.tol=.001)))

  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., ind)) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#### Hierarchical: Function to combine both model fit stats for all levels into one dataframe ####

multi_df_hier <- function(model,n,estimator,reps){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- multi_fit_hier(model,n,estimator,reps)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_hier(model,n,estimator,reps)

  #Produce final table of fit indices for each level (as a list)
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}

############################################################################
################################# hier2 ####################################
############################################################################

#### Function for Isolating Second-Order Factors ####

iso_second <- function(model){

  #isolate first-order factors
  factFirst <- iso_first(model)

  #isolate second-order factors
  factSecond <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(op=="=~") %>%
    dplyr::filter(rhs %in% factFirst) %>%
    dplyr::select(lhs) %>%
    base::unique() %>%
    base::unlist()

  return(factSecond)
}

#### Function to Create lav_file with Only Second-Order Factors ####

lav_file_second <- function(model){

  #isolate names of first-order factors
  factFirst <- iso_first(model)

  #isolate names of second-order factors
  factSecond <- iso_second(model)

  #separate out factor correlation rows
  corr_lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(op=="~~") %>%
    dplyr::filter(rhs != lhs) %>%
    dplyr::filter(lhs %in% factSecond) %>%
    dplyr::filter(rhs %in% factSecond)

  #isolate second-order factor rows
  fact_lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(op=="=~") %>%
    dplyr::filter(lhs %in% factSecond) %>%
    dplyr::filter(rhs %in% factFirst)

  #prep the model
  lav_file <- dplyr::bind_rows(fact_lav_file, corr_lav_file)

  return(lav_file)
}

### Hierarchical2: Function to see which items are available ###

one_num_hier2 <- function(model){

  #Rename (just to be consistent with shiny app)
  Mod_C <- model

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lav_file_second(Mod_C)

  #identify all factor names
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #Identify any items that already have an error covariance
  items_covariance <- factors %>%
    dplyr::mutate(type="Factor") %>%
    dplyr::full_join(lav_file, by = "lhs") %>%
    dplyr::select(-type,type) %>%
    dplyr::select(lhs,op,rhs,type) %>%
    dplyr::filter(op=="=~" | is.na(type)) %>%
    dplyr::filter(is.na(type)) %>%
    dplyr::select(-type) %>%
    tidyr::pivot_longer(-op,names_to = "test", values_to = "rhs") %>%
    dplyr::select(-op,-test) %>%
    dplyr::mutate(lhs=NA,op=NA,ustart=NA)

  #Isolate the items that do not already have an error covariance
  solo_items <- lav_file %>%
    dplyr::select(lhs,op,rhs,ustart) %>%
    base::rbind(items_covariance) %>%
    dplyr::filter(op=="=~"|is.na(op)) %>%
    dplyr::group_by(rhs) %>%
    dplyr::add_tally() %>%
    dplyr::filter(n==1) %>%
    dplyr::ungroup() %>%
    arrange(abs(ustart))

  return(solo_items)
}

#### Hierarchical2: Function to create misspecification statement ####

one_add_hier2 <- function(model){

  #Read in available items
  itemoptions <- one_num_hier2(model)

  #Count number of available items
  num_i <- base::nrow(itemoptions)

  #Select items for misspecification depending on number of available items
  if(num_i==4){
    num_m <- itemoptions %>%
      dplyr::slice(1:2)
  }else if(num_i==5){
    num_m <- itemoptions %>%
      dplyr::slice(1:4)
  }else{
    num_m <- itemoptions %>%
      dplyr::slice(1:(floor(num_i/2)*2))
  }

  #Identifiers to separate odds and even rows
  evenindex <- base::seq(2,base::nrow(num_m),2)
  oddindex <- base::seq(1,base::nrow(num_m),2)

  #Separate
  left <- num_m[evenindex,]
  right <- num_m[oddindex,] %>%
    `colnames<-`(c("lhs_1","op_1","rhs_1","ustart_1","n_1"))

  #Create misspecification statements
  Residual_Correlation <- base::cbind(left,right) %>%
    dplyr::mutate(cor=.5,
                  opp="~~",
                  star="*") %>%
    tidyr::unite(V1,c("rhs","opp","cor","star","rhs_1"),sep=" ") %>%
    dplyr::select(V1)

  return(Residual_Correlation)
}


#### Hierarchical2: Function to create Misspecified DGM ####

DGM_one_hier2 <- function(model){

  #Count number of available items for number of misspecifications
  num_m<- base::nrow(one_num_hier2(model))

  #Figure out number of levels given number of available items
  if(num_m==4){
    L1 <- 1
    L2 <- 2
    levels <- base::rbind(L1,L2)
  }else if(num_m==5){
    L1 <- 1
    L2 <- 2
    levels <- base::rbind(L1,L2)
  }else{
    L3 <- base::floor(num_m/2)
    L2 <- base::floor((2*L3)/3)
    L1 <- base::floor(L3/3)
    levels <- base::rbind(L1,L2,L3)
  }

  #Read in misspecifications
  mod <- one_add_hier2(model)

  #Get parameters for true dgm
  Mods <- model
  #Mod_C <- base::as.character(Mods$V1)

  #single_mod <- base::lapply(levels, function(x) base::rbind(Mod_C,mod[base::seq(x), ,drop = FALSE]) %>%
  #                       base::data.frame() %>%
  #                       dplyr::pull(V1))
  #This made you miserable. Shiny was struggling with \n at the end of strings here, for some reason.

  #Create a list for every row in the mod object (misspecifications)
  #For each element, bind the misspecification to the OG model statement sequentially
  #Turn it into a dataframe and extract
  single_mod <- base::lapply(levels, function(x) base::rbind(Mods,mod[base::seq(x), ,drop = FALSE]) %>%
                               base::data.frame() %>%
                               dplyr::pull(V1))

  return(single_mod)

}

### Hierarchical2: Simulate fit indices for misspecified model for all levels ###

one_fit_hier2 <- function(model,n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm
  misspec_dgm <- DGM_one_hier2(model)

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set seed
  set.seed(649364)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Simulate one large dataset for each misspecification (use map to apply across each
  #element (set of misspecifications) in the list)
  all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*r,
                                                                            latent=FALSE,errors=FALSE))

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_misspec <- purrr::map(all_data_misspec,~base::cbind(.,rep_id_misspec))

  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~dplyr::group_by(.,rep_id_misspec) %>%
                               tidyr::nest())

  #Grab data level of the list
  data <- purrr::map(misspec_data,2)

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa for each element in the list
  misspec_cfa <- purrr::map(data, function(x) purrr::map(x, function(y) lavaan::cfa(model=mod,
                                                                                    estimator=estimator,
                                                                                    data=y,
                                                                                    std.lv=TRUE,
                                                                                    se=se,
                                                                                    check.gradient=FALSE,
                                                                                    check.post=FALSE,
                                                                                    check.vcov=FALSE,
                                                                                    control=list(rel.tol=.001))))

  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, ind)) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))

  set.seed(NULL)

  return(misspec_fit_sum)

}

#### Hierarchical2: Function to create True DGM (aka, just the model the user read in) ####

true_fit_one_hier2 <- function(model,n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  true_dgm <- model

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set Seed
  set.seed(326267)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Simulate one large dataset
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n*r,
                                                 latent = FALSE,
                                                 errors = FALSE)

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_true <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)

  #Group and list
  true_data <- dat_rep_true %>%
    dplyr::group_by(rep_id_true) %>%
    tidyr::nest() %>%
    base::as.list()

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa
  true_cfa <- purrr::map(true_data$data,function(x) lavaan::cfa(model=mod,
                                                                estimator=estimator,
                                                                data=x,
                                                                std.lv=TRUE,
                                                                se=se,
                                                                check.gradient=FALSE,
                                                                check.post=FALSE,
                                                                check.vcov=FALSE,
                                                                control=list(rel.tol=.001)))

  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., ind)) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#### Hierarchical2: Function to combine both model fit stats for all levels into one dataframe ####

one_df_hier2 <- function(model,n,estimator,reps){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- one_fit_hier2(model,n,estimator,reps)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_one_hier2(model,n,estimator,reps)

  #Produce final table by level
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}

##################################################
################# catOne #########################
##################################################

###### Separate out thresholds from manual=TRUE#
##### this is used to to create threshold list for simulation discretized #
cleanthreshold <- function(model){

  suppressMessages(model %>%
                     lavaan::lavaanify(fixed.x = FALSE) %>%
                     dplyr::filter(grepl("^t",rhs)) %>%
                     dplyr::select(lhs,op,rhs,ustart) %>%
                     dplyr::rename(est.std = ustart))
}

### discard thresholds if Manual = TRUE but keep the estimates#
### Basiclly a long way to get equivalent of model statement with manual=TRUE for continuous case#
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

### function to grab the thresholds if input is a lavaan file
Thresh <- function(model){

  #Extract standardized solution from lavaan object
  lav <- lavaan::standardizedSolution(model)

  #Create model statement
  t <- suppressMessages(as.data.frame(lav %>%
                                        dplyr::filter(lhs != rhs) %>%
                                        dplyr::group_by(lhs,op) %>%
                                        dplyr::filter(op == "|") %>%
                                        dplyr::select(lhs,op,rhs,est.std)))

  return(t)

}


##function to identify items with thresholds
##Used for identifying which items are categorical in simulation
##That allows for a mix of continuous and categorical items, if someone happens to have that
cat_items <- function(threshold){
  #prep the model
  items <- threshold %>%
    dplyr::filter(op=="|") %>%
    dplyr::select(lhs) %>%
    base::unique()

  return(items)
}

### return thresholds for all categorical items##
## works for lavaan input or manual = TRUE ##
## cleanthreshold function puts thresholds into same format as lavaan object ##
## items with different number of categories result in NA for extra categories##
th<-function(threshold){
  a1<-stats::reshape(threshold, idvar="rhs",timevar="lhs",direction="wide", drop="op") #transpose so that each column is var and threshold read down
  a1<-a1[-c(1)] #remove unnecessary columns
  a2<-t(as.data.frame(cat_items(threshold))) #transpose names to row vector
  names(a1)<-c(a2)#replace names so that they align
  a2<-as.data.frame(a2)
  return(a1)
}
### return item names for categorical items, this allows continuous & categorical in same model ###
th2<-function(threshold){
  a1<-stats::reshape(threshold, idvar="rhs",timevar="lhs",direction="wide", drop="op") #transpose so that each column is var and threshold read down
  a1<-a1[-c(1)] #remove unnecessary columns
  a2<-t(as.data.frame(cat_items(threshold))) #transpose names to row vector
  names(a1)<-c(a2)#replace names so that they align
  a2<-as.data.frame(a2)
  return(a2)
}

one_fit_cat <- function(model,n,reps,threshold, estimator){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm
  misspec_dgm <- DGM_one(model)

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set seed
  set.seed(649364)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Simulate one large dataset for each misspecification (use map to apply across each
  #element (set of misspecifications) in the list)
  all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*r,
                                                                            latent=FALSE,errors=FALSE))
  #Loop through each categorical items in the simulated data
  # recode simulated MVN data into categories
  #multiply by 100 to space categories out and avoid overwriting (category labels are ordinal, exact label irrelevant)
  ### extra "x" loop here to loop through each Level of misspecification
  ### "x" loop not present in true data generation function
  a1<-th(threshold)
  a2<-th2(threshold)
  for (x in 1:length(misspec_dgm))
  {
    for (i in 1:ncol(a1))
    {
      u<- as.numeric(sum(!is.na(a1[,i])))
      all_data_misspec[[x]]<-all_data_misspec[[x]] %>%
        dplyr::mutate(!!a2[,i] := case_when(
          !!rlang::sym(a2[,i])  <= a1[1,i] ~100,
          !!rlang::sym(a2[,i]) >   a1[u,i] ~100*(u+1),
          TRUE ~ !!rlang::sym(a2[,i]))
        )
      if(sum(!is.na(a1[,i])) > 1){
        for (j in 1:sum(!is.na(a1[,i]))) {
          all_data_misspec[[x]]<-all_data_misspec[[x]] %>%
            dplyr::mutate(!!a2[,i] := case_when(
              between( !!rlang::sym(a2[,i]) , a1[j,i], a1[j+1,i]) ~ as.numeric(100*(j+1)),
              TRUE~ !!rlang::sym(a2[,i]))
            )
        }
      }
    }
  }

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_misspec <- purrr::map(all_data_misspec,~base::cbind(.,rep_id_misspec))

  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~dplyr::group_by(.,rep_id_misspec) %>%
                               tidyr::nest())

  #Grab data level of the list
  data <- purrr::map(misspec_data,2)

  #categorical vector of categorical items
  clist<-cat_items(threshold)[,1]

  se<-ifelse(startsWith(estimator,"ULS"),"robust.sem","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa
  #changed lavaan code to treat items in clist as categorical
  #also added some options to suppress common warnings
  misspec_cfa <- purrr::map(data, function(x) purrr::map(x, function(y) lavaan::cfa(model = mod, data=y, estimator=estimator, ordered=clist, std.lv=TRUE, se=se,
                                                                                    check.gradient=FALSE,
                                                                                    check.post=FALSE,
                                                                                    check.vcov=FALSE,
                                                                                    control=list(rel.tol=.001))))


  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y,ind)) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))

  set.seed(NULL)

  return(misspec_fit_sum)
  }

#simulates MVN data, but then bins it based on thresholds
true_fit_one_cat <- function(model,n,reps, threshold, estimator){

  #Get clean model equation
  mod <- cleanmodel(model)

  true_dgm <- model

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set Seed
  set.seed(326267)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Simulate one large dataset
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n*r,
                                                 latent = FALSE,
                                                 errors = FALSE)
  #similar to above, loops through categorical items to discretize data based on thresholds
  a1<-th(threshold)
  a2<-th2(threshold)

  for (i in 1:ncol(a1))
  {
    u<- as.numeric(sum(!is.na(a1[,i])))
    all_data_true<-all_data_true%>%
      dplyr::mutate(!!a2[,i] := case_when(
        !!rlang::sym(a2[,i])  <= a1[1,i] ~100,
        !!rlang::sym(a2[,i]) >   a1[u,i] ~100*(u+1),
        TRUE ~ !!rlang::sym(a2[,i]))
      )
    if(sum(!is.na(a1[,i])) > 1){
      for (j in 1:sum(!is.na(a1[,i]))) {
        all_data_true<-all_data_true%>%
          dplyr::mutate(!!a2[,i] := case_when(
            between( !!rlang::sym(a2[,i]) , a1[j,i], a1[j+1,i]) ~ as.numeric(100*(j+1)),
            TRUE~ !!rlang::sym(a2[,i]))
          )
      }
    }
  }


  #Create indicator to split into 500 datasets for 500 reps
  rep_id_true <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)

  #Group and list
  true_data <- dat_rep_true %>%
    dplyr::group_by(rep_id_true) %>%
    tidyr::nest() %>%
    base::as.list()

   #create character vector listing categorical items
  clist<-cat_items(threshold)[,1]

  se<-ifelse(startsWith(estimator,"ULS"),"robust.sem","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa
  #added ORDERED = to treat items as categorical
  #also added options to suppress common warnings
  true_cfa <- purrr::map(true_data$data,function(x) lavaan::cfa(model = mod, data=x, std.lv=TRUE, ordered=clist, estimator=estimator, std.lv=TRUE, se=se,
                                                                check.gradient=FALSE,
                                                                check.post=FALSE,
                                                                check.vcov=FALSE,control=list(rel.tol=.001)))

  #Extract fit stats from each rep (list) into a data frame and clean
  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., ind)) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)
  }

one_df_cat <- function(model,n,reps,threshold, estimator){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- one_fit_cat(model,n,reps, threshold,estimator)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_one_cat(model,n,reps, threshold,estimator)

  #Produce final table by level
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}

##################################################
################# catHB ##########################
##################################################

multi_fit_HB_cat <- function(model,n,reps, threshold, estimator){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm (this is a list)
  misspec_dgm <- DGM_Multi_HB(model)

  #Use max sample size of 2000
  n <- min(n,2000)

  #Set seed
  set.seed(269854)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Simulate one large dataset for each misspecification (use map to apply across each
  #element (set of misspecifications) in the list)
  all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*r,latent=FALSE,errors=FALSE))

  a1<-th(threshold)
  a2<-th2(threshold)
  for (x in 1:length(misspec_dgm))
  {
    for (i in 1:ncol(a1))
    {
      u<- as.numeric(sum(!is.na(a1[,i])))
      all_data_misspec[[x]]<-all_data_misspec[[x]] %>%
        dplyr::mutate(!!a2[,i] := case_when(
          !!rlang::sym(a2[,i])  <= a1[1,i] ~100,
          !!rlang::sym(a2[,i]) >   a1[u,i] ~100*(u+1),
          TRUE ~ !!rlang::sym(a2[,i]))
        )
      if(sum(!is.na(a1[,i])) > 1){
        for (j in 1:sum(!is.na(a1[,i]))) {
          all_data_misspec[[x]]<-all_data_misspec[[x]] %>%
            dplyr::mutate(!!a2[,i] := case_when(
              between( !!rlang::sym(a2[,i]) , a1[j,i], a1[j+1,i]) ~ as.numeric(100*(j+1)),
              TRUE~ !!rlang::sym(a2[,i]))
            )
        }
      }
    }
  }

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- rep(1:r,n)

  #Combine indicator with dataset for each element in list
  dat_rep_misspec <- purrr::map(all_data_misspec,~cbind(.,rep_id_misspec))

  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~group_by(.,rep_id_misspec) %>%
                               tidyr::nest())

  #Grab data level of the list
  data <- purrr::map(misspec_data,2)

  #categorical vector of categorical items
  clist<-cat_items(threshold)[,1]

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"robust.sem","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa for each element in the list
  misspec_cfa <- purrr::map(data, function(x) purrr::map(x, function(y) lavaan::cfa(model=mod,
                                                                                    estimator=estimator,
                                                                                    data=y,
                                                                                    ordered=clist,
                                                                                    std.lv=TRUE,
                                                                                    se=se,
                                                                                    check.gradient=FALSE,
                                                                                    check.post=FALSE,
                                                                                    check.vcov=FALSE,
                                                                                    control=list(rel.tol=.001))))

  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, ind)) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))
  set.seed(NULL)

  return(misspec_fit_sum)

}

true_fit_HB_cat <- function(model,n,reps, threshold, estimator){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for true DGM
  true_dgm <- model

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Set Seed
  set.seed(267326)

  #Simulate one large dataset
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n*r,
                                                 latent = FALSE,
                                                 errors = FALSE)

  a1<-th(threshold)
  a2<-th2(threshold)

  for (i in 1:ncol(a1))
  {
    u<- as.numeric(sum(!is.na(a1[,i])))
    all_data_true<-all_data_true%>%
      dplyr::mutate(!!a2[,i] := case_when(
        !!rlang::sym(a2[,i])  <= a1[1,i] ~100,
        !!rlang::sym(a2[,i]) >   a1[u,i] ~100*(u+1),
        TRUE ~ !!rlang::sym(a2[,i]))
      )
    if(sum(!is.na(a1[,i])) > 1){
      for (j in 1:sum(!is.na(a1[,i]))) {
        all_data_true<-all_data_true%>%
          dplyr::mutate(!!a2[,i] := case_when(
            between( !!rlang::sym(a2[,i]) , a1[j,i], a1[j+1,i]) ~ as.numeric(100*(j+1)),
            TRUE~ !!rlang::sym(a2[,i]))
          )
      }
    }
  }

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_true <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)

  #Group and list
  true_data <- dat_rep_true %>%
    dplyr::group_by(rep_id_true) %>%
    tidyr::nest() %>%
    base::as.list()

  #create character vector listing categorical items
  clist<-cat_items(threshold)[,1]

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"robust.sem","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa
  true_cfa <- purrr::map(true_data$data,function(x) lavaan::cfa(model=mod,
                                                                estimator=estimator,
                                                                data=x,
                                                                ordered=clist,
                                                                std.lv=TRUE,
                                                                se=se,
                                                                check.gradient=FALSE,
                                                                check.post=FALSE,
                                                                check.vcov=FALSE,
                                                                control=list(rel.tol=.001)))

  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., ind)) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#### Multi-Factor: Function to combine both model fit stats for all levels into one dataframe ####

multi_df_HB_cat <- function(model,n,reps, threshold,estimator){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- multi_fit_HB_cat(model,n,reps,threshold, estimator)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_HB_cat(model,n,reps,threshold, estimator)

  #Produce final table of fit indices for each level (as a list)
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}

##################################################
################# nnorOne ########################
##################################################

### One-factor: Simulate fit indices for misspecified model for all levels ###
one_fit_nnor <- function(model,data, n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm
  misspec_dgm <- DGM_one(model)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  miss<-list()
  for(i in 1:length(misspec_dgm))
  {
    dat <- simstandard::sim_standardized_matrices(misspec_dgm[[i]])
    x<-dat$Correlations$R
    l<-nrow(x)
    a<-x[1:(l-1),1:(l-1)]
    a<-a[order(rownames(a)),order(colnames(a))]
    ll<-nrow(a)
    names<-colnames(a)

    miss[[i]]<-a
  }

  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]

  n <- base::min(nrow(data1),2000)

  #vector of all 0s
  mu<-c(colMeans(data1,na.rm=T))



  data_m<-list()
  for(i in 1:length(misspec_dgm))
  {
    data_m[[i]]<-semTools::bsBootMiss(Sigma =miss[[i]], Mu=mu, rawData = data1, nBoot=reps,bootSamplesOnly  = TRUE, seed=649364)
  }

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa for each element in the list
  misspec_cfa <- purrr::map(data_m, function(x) purrr::map(x, function(y) lavaan::cfa(model=mod,
                                                                                      estimator=estimator,
                                                                                      data=y,
                                                                                      std.lv=TRUE,
                                                                                      se=se,
                                                                                      check.gradient=FALSE,
                                                                                      check.post=FALSE,
                                                                                      check.vcov=FALSE,
                                                                                      control=list(rel.tol=.001))))

  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, ind)) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))

  set.seed(NULL)

  return(misspec_fit_sum)

}

#### One_Factor: Function to create True DGM (aka, just the model the user read in) ####

true_fit_one_nnor <- function(model,data,n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  true_dgm <- model

  #Set Seed
  set.seed(326267)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  nf<-number_factor(model)
  dat <- simstandard::sim_standardized_matrices(true_dgm)
  x<-dat$Correlations$R
  l<-nrow(x)
  a<-x[1:(l-nf),1:(l-nf)]
  a<-a[order(rownames(a)),order(colnames(a))]
  names<-colnames(a)
  true<-a


  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]

  n <- base::min(nrow(data1),2000)

  #vector of all 0s
  mu<-c(colMeans(data1,na.rm=T))


  data_t<-semTools::bsBootMiss(Sigma =true, Mu=mu, rawData = data1, nBoot=reps,bootSamplesOnly  = TRUE)


  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa
  true_cfa <- purrr::map(data_t,function(x) lavaan::cfa(model=mod,
                                                        estimator=estimator,
                                                        data=x,
                                                        std.lv=TRUE,
                                                        se=se,
                                                        check.gradient=FALSE,
                                                        check.post=FALSE,
                                                        check.vcov=FALSE,
                                                        control=list(rel.tol=.001)))

  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., ind)) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#### One-Factor: Function to combine both model fit stats for all levels into one dataframe ####

one_df_nnor <- function(model,data,n,estimator,reps){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- one_fit_nnor(model,data,n,estimator,reps)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_one_nnor(model,data,n,estimator,reps)

  #Produce final table by level
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}

data_nnor <- function(model,data, n,reps){

  #Get parameters for misspecified dgm
  mod <- cleanmodel(model)

  true_dgm <- model

  #Set Seed
  set.seed(326267)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  nf<-number_factor(model)
  dat <- simstandard::sim_standardized_matrices(true_dgm)
  x<-dat$Correlations$R
  l<-nrow(x)
  a<-x[1:(l-nf),1:(l-nf)]
  a<-a[order(rownames(a)),order(colnames(a))]
  names<-colnames(a)
  true<-a


  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]

  n <- base::min(nrow(data1),2000)

  #vector of all 0s
  mu<-c(colMeans(data1,na.rm=T))

  data_t<-semTools::bsBootMiss(Sigma =true, Mu=mu, rawData = data1, nBoot=reps,bootSamplesOnly  = TRUE)

  ##Aggregate all bootstrapped data into one matrix per level
  data_t_flat<-list()

  data_t_flat<-base::do.call(rbind, data_t)

  miss_dat<-list()

  miss_dat$data_t<-data_t_flat
  miss_dat$data1<-data1

  return(miss_dat)
}

##################################################
################# nnorHB #########################
##################################################

### multi-factor: Simulate fit indices for misspecified model for all levels ###
multi_fit_nnor <- function(model,data, n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm
  misspec_dgm <- DGM_Multi_HB(model)

  #Number of reps
  r <- reps

  #numbber of factors
  nf<-number_factor(model)

  #generate model-implied matrices for misspecified models
  miss<-list()
  for(i in 1:length(misspec_dgm))
  {
    dat <- simstandard::sim_standardized_matrices(misspec_dgm[[i]])
    x<-dat$Correlations$R
    l<-nrow(x)
    a<-x[1:(l-nf),1:(l-nf)]
    a<-a[order(rownames(a)),order(colnames(a))]
    ll<-nrow(a)
    names<-colnames(a)

    miss[[i]]<-a
  }

  #identify names of variables from dataset
  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]

  n <- base::min(nrow(data1),2000)

  #vector of all 0s
  mu<-c(colMeans(data1,na.rm=T))

  #Bollen-Stein Bootstrap data based on misspecified model-implied matrices
  data_m<-list()
  for(i in 1:length(misspec_dgm))
  {
    data_m[[i]]<-semTools::bsBootMiss(Sigma =miss[[i]], Mu=mu, rawData = data1, nBoot=reps,bootSamplesOnly  = TRUE, seed=649364)
  }

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa for each element in the list
  misspec_cfa <- purrr::map(data_m, function(x) purrr::map(x, function(y) lavaan::cfa(model=mod,
                                                                                      estimator=estimator,
                                                                                      data=y,
                                                                                      std.lv=TRUE,
                                                                                      se=se,
                                                                                      check.gradient=FALSE,
                                                                                      check.post=FALSE,
                                                                                      check.vcov=FALSE,
                                                                                      control=list(rel.tol=.001))))

  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, ind)) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))

  set.seed(NULL)

  return(misspec_fit_sum)

}

#### multi_Factor: Function to create True DGM (aka, just the model the user read in) ####

true_fit_multi_nnor <- function(model,data,n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  true_dgm <- model

  #Set Seed
  set.seed(326267)

  #Number of reps
  r <- reps

  #number of factors
  nf<-number_factor(model)

  #model implied correlation for model consistent with user's model
  dat <- simstandard::sim_standardized_matrices(true_dgm)
  x<-dat$Correlations$R
  l<-nrow(x)
  a<-x[1:(l-nf),1:(l-nf)]
  a<-a[order(rownames(a)),order(colnames(a))]
  names<-colnames(a)
  true<-a

  #select variables from the data that were used in the model
  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]

  n <- base::min(nrow(data1),2000)

  #vector of all 0s
  mu<-c(colMeans(data1,na.rm=T))

  #bootstrap data that are consistent with the user's model
  data_t<-semTools::bsBootMiss(Sigma =true, Mu=mu, rawData = data1, nBoot=reps,bootSamplesOnly  = TRUE)

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa
  true_cfa <- purrr::map(data_t,function(x) lavaan::cfa(model=mod,
                                                        estimator=estimator,
                                                        data=x,
                                                        std.lv=TRUE,
                                                        se=se,
                                                        check.gradient=FALSE,
                                                        check.post=FALSE,
                                                        check.vcov=FALSE,
                                                        control=list(rel.tol=.001)))

  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., ind)) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#### multi-Factor: Function to combine both model fit stats for all levels into one dataframe ####

multi_df_nnor <- function(model,data,n,estimator,reps){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- multi_fit_nnor(model,data,n,estimator,reps)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_multi_nnor(model,data,n,estimator,reps)

  #Produce final table by level
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}


##################################################
################# likertHB #######################
##################################################


multi_fit_likert <- function(model,data, n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm
  misspec_dgm <- DGM_Multi_HB(model)

  n <- min(n,2000)

  #Set seed
  set.seed(269854)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #numbber of factors
  nf<-number_factor(model)

  #Simulate one large dataset for each misspecification (use map to apply across each
  #element (set of misspecifications) in the list)
  all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*r,latent=FALSE,errors=FALSE))

  names<-colnames(all_data_misspec[[1]])

  #identify names of variables from dataset
  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]


  #create empty matrix for proportions (g) and threshold (p)
  g<-matrix(nrow=(max(data1,na.rm=T)-min(data1, na.rm=T)), ncol=ncol(data1))
  p<-g

  #loop through all variables in fitted model
  for (i in 1:ncol(data1)){
    #loop through all categories in fitted model
    for (h in min(data1[,i], na.rm=T):(max(data1[,i],na.rm=T)-1)){
      #proportion of responses at or below each category
      g[h,i]<-sum(table(data1[,i])[names(table(data1[,i]))<=h])/nrow(data1)
      #inverse standard normal to determine thresholds
      p[h,i]<-qnorm(c(g[h,i]))
      dp<-as.data.frame(p)
    }
  }

  colnames(dp)<-colnames(data1)
  a2<-colnames(dp)

  for (i in 1:ncol(data1)){
    if(is.na(dp[1,i])){
      dp[1:min(which(!is.na(dp[,i]))-1),i]=-10
    }
  }


  #loop through misspecifications
  for (x in 1:length(misspec_dgm))
  {
    for (i in 1:ncol(dp))
    {# first loop transforms highest category
      u<- as.numeric(sum(!is.na(dp[,i])))
      all_data_misspec[[x]]<-all_data_misspec[[x]] %>%
        dplyr::mutate(!!a2[i] := case_when(
          !!rlang::sym(a2[i])  <= dp[1,i] ~100,
          !!rlang::sym(a2[i]) >   dp[u,i] ~100*(u+1),
          TRUE ~ !!rlang::sym(a2[i]))
        )
      #second loop transforms all lower categories
      if(sum(!is.na(dp[,i])) > 1){
        for (j in 1:sum(!is.na(dp[,i]))) {
          all_data_misspec[[x]]<-all_data_misspec[[x]] %>%
            dplyr::mutate(!!a2[i] := case_when(
              between( !!rlang::sym(a2[i]) , dp[j,i], dp[j+1,i]) ~ as.numeric(100*(j+1)),
              TRUE~ !!rlang::sym(a2[i]))
            )
        }
      }
    }
    #loops multiply by 100 to avoid overwriting MVN data with values above 1
    #divide by 100 to put things back onto Likert metric
    all_data_misspec[[x]]<-all_data_misspec[[x]]/100
  }


  rep_id_misspec <- rep(1:r,n)

  #Combine indicator with dataset for each element in list
  dat_rep_misspec <- purrr::map(all_data_misspec,~cbind(.,rep_id_misspec))

  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~group_by(.,rep_id_misspec) %>%
                               tidyr::nest())

  #Grab data level of the list
  data_m <- purrr::map(misspec_data,2)

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa for each element in the list
  misspec_cfa <- purrr::map(data_m, function(x) purrr::map(x, function(y) lavaan::cfa(model=mod,
                                                                                      estimator=estimator,
                                                                                      data=y,
                                                                                      std.lv=TRUE,
                                                                                      se=se,
                                                                                      check.gradient=FALSE,
                                                                                      check.post=FALSE,
                                                                                      check.vcov=FALSE,
                                                                                      control=list(rel.tol=.001))))

  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, ind)) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))

  set.seed(NULL)

  return(misspec_fit_sum)

}

#### multi_Factor: Function to create True DGM (aka, just the model the user read in) ####

true_fit_multi_likert <- function(model,data,n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  true_dgm <- model

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set Seed
  set.seed(326267)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Simulate one large dataset
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n*r,
                                                 latent = FALSE,
                                                 errors = FALSE)

  names<-colnames(all_data_true)

  #identify names of variables from dataset
  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]

  #create empty matrix for proportions (g) and threshold (p)
  g<-matrix(nrow=(max(data1, na.rm=T)-min(data1, na.rm=T)), ncol=ncol(data1))
  p<-g

  #loop through all variables in fitted model
  for (i in 1:ncol(data1)){
    #loop through all categories in fitted model
    for (h in min(data1[,i], na.rm=T):(max(data1[,i], na.rm=T)-1)){
      #proportion of responses at or below each category
      g[h,i]<-sum(table(data1[,i])[names(table(data1[,i]))<=h])/nrow(data1)
      #inverse standard normal to determine thresholds
      p[h,i]<-qnorm(c(g[h,i]))
      dp<-as.data.frame(p)
    }
  }

  colnames(dp)<-colnames(data1)
  a2<-colnames(dp)

  for (i in 1:ncol(data1)){
  if(is.na(dp[1,i])){
  dp[1:min(which(!is.na(dp[,i]))-1),i]=-10
  }
  }

      for (i in 1:ncol(dp))
  {# first loop transforms highest category
    u<- as.numeric(sum(!is.na(dp[,i])))
    all_data_true<-all_data_true %>%
      dplyr::mutate(!!a2[i] := case_when(
        !!rlang::sym(a2[i])  <= dp[1,i] ~100,
        !!rlang::sym(a2[i]) >   dp[u,i] ~100*(u+1),
        TRUE ~ !!rlang::sym(a2[i]))
      )
    #second loop transforms all lower categories
    if(sum(!is.na(dp[,i])) > 1){
      for (j in 1:sum(!is.na(dp[,i]))) {
        all_data_true<-all_data_true %>%
          dplyr::mutate(!!a2[i] := case_when(
            between( !!rlang::sym(a2[i]) , dp[j,i], dp[j+1,i]) ~ as.numeric(100*(j+1)),
            TRUE~ !!rlang::sym(a2[i]))
          )
      }
    }
  }
  #loops multiply by 100 to avoid overwriting MVN data with values above 1
  #divide by 100 to put things back onto Likert metric
  all_data_true<-all_data_true/100

  rep_id_true <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)

  #Group and list
  data_t <- dat_rep_true %>%
    dplyr::group_by(rep_id_true) %>%
    tidyr::nest() %>%
    base::as.list()


  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa
  true_cfa <- purrr::map(data_t$data,function(x) lavaan::cfa(model=mod,
                                                             estimator=estimator,
                                                             data=x,
                                                             std.lv=TRUE,
                                                             se=se,
                                                             check.gradient=FALSE,
                                                             check.post=FALSE,
                                                             check.vcov=FALSE,
                                                             control=list(rel.tol=.001)))

  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., ind)) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#### multi-Factor: Function to combine both model fit stats for all levels into one dataframe ####

multi_df_likert <- function(model,data,n,estimator,reps){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- multi_fit_likert(model,data,n,estimator,reps)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_multi_likert(model,data,n,estimator,reps)

  #Produce final table by level
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}

data_likert <- function(model,data, n){

  #Get clean model equation
  mod <- cleanmodel(model)

  true_dgm <- model

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set Seed
  set.seed(326267)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)

  #Simulate one large dataset
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n,
                                                 latent = FALSE,
                                                 errors = FALSE)

  names<-colnames(all_data_true)

  #identify names of variables from dataset
  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]

  #create empty matrix for proportions (g) and threshold (p)
  g<-matrix(nrow=(max(data1, na.rm=T)-min(data1, na.rm=T)), ncol=ncol(data1))
  p<-g

  #loop through all variables in fitted model
  for (i in 1:ncol(data1)){
    #loop through all categories in fitted model
    for (h in min(data1[,i], na.rm=T):(max(data1[,i], na.rm=T)-1)){
      #proportion of responses at or below each category
      g[h,i]<-sum(table(data1[,i])[names(table(data1[,i]))<=h])/nrow(data1)
      #inverse standard normal to determine thresholds
      p[h,i]<-qnorm(c(g[h,i]))
      dp<-as.data.frame(p)
    }
  }

  colnames(dp)<-colnames(data1)
  a2<-colnames(dp)

  for (i in 1:ncol(data1)){
    if(is.na(dp[1,i])){
      dp[1:min(which(!is.na(dp[,i]))-1),i]=-10
    }
  }

  for (i in 1:ncol(dp))
  {# first loop transforms highest category
    u<- as.numeric(sum(!is.na(dp[,i])))
    all_data_true<-all_data_true %>%
      dplyr::mutate(!!a2[i] := case_when(
        !!rlang::sym(a2[i])  <= dp[1,i] ~100,
        !!rlang::sym(a2[i]) >   dp[u,i] ~100*(u+1),
        TRUE ~ !!rlang::sym(a2[i]))
      )
    #second loop transforms all lower categories
    if(sum(!is.na(dp[,i])) > 1){
      for (j in 1:sum(!is.na(dp[,i]))) {
        all_data_true<-all_data_true %>%
          dplyr::mutate(!!a2[i] := case_when(
            between( !!rlang::sym(a2[i]) , dp[j,i], dp[j+1,i]) ~ as.numeric(100*(j+1)),
            TRUE~ !!rlang::sym(a2[i]))
          )
      }
    }
  }
  #loops multiply by 100 to avoid overwriting MVN data with values above 1
  #divide by 100 to put things back onto Likert metric
  all_data_true<-all_data_true/100

  miss_dat<-list()

  miss_dat$sim<-all_data_true
  miss_dat$orig<-data1

  return(miss_dat)
}

##################################################
################# likertHB2 ######################
##################################################


#function for misspecified model
multi_fit_likert2 <- function(model,data, n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  factors <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(lhs != rhs)%>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #Get model-implied matrix of input model
  norm <- simstandard::sim_standardized_matrices(model)
  normx<-norm$Correlations$R
  l<-nrow(normx)-nrow(factors)
  a<-normx[1:l,1:l]
  diag(a)=1

  #names of variable in model
  names<-colnames(a)

  #identify names of variables from dataset
  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]

  #rescale so that minimum value is always ==1
  #needed to properly index computations below
  #also needed to avoid 0s because there is multiplication involved
  #will be scaled back at the end after loop indexes are not longer needed
  d2<-matrix(sapply(data1, function(x) min(x, na.rm=T)-1), nrow=1, ncol=ncol(data1))
  d3<-matrix(rep(d2,each=nrow(data1)), nrow=nrow(data1), ncol=ncol(data1))
  d4<-matrix(rep(d2,each=nrow(data1)), nrow=r*nrow(data1), ncol=ncol(data1))
  colnames(d3)<-colnames(data1)
  colnames(d4)<-colnames(data1)
  data1<-data1-d3

  #create empty matrix for proportions in each category (g)
  #GenOrd requires list (p) for each item
  #t is thresholds for eventual discretization
  g<-matrix(nrow=(max(data1,na.rm=T)-min(data1, na.rm=T)), ncol=ncol(data1))
  p<-list()
  t<-g

  #loop through all variables in fitted model
  for (i in 1:ncol(data1)){
    #loop through all categories for each variable
    for (h in min(data1[,i],na.rm=T):(max(data1[,i],na.rm=T)-1)){
      #proportion of responses at or below each category, only count non-missing in the denominator
      g[h,i]<-sum(table(data1[,i])[names(table(data1[,i]))<=h])/sum(!is.na(data1[,i]))
      p[[i]]<-g[,i]
      t[h,i]<-qnorm(c(g[h,i]))
      dp<-as.data.frame(t)
    }
  }

  #handle different number of categories per item
  p1<-list()
  for (x in 1:length(p)){
    p1[[x]]<-p[[x]][!is.na(p[[x]])]
  }

  #Assign column names so it's clear which thresholds go with which variable
  colnames(g)<-colnames(data1)
  colnames(dp)<-colnames(data1)
  a2<-colnames(g)

  #misspecifications on observed scale
  misspec_dgm <- DGM_Multi_HB(model)

  #model-implied correlation after discretization
  Miss_Cor <- purrr::map(misspec_dgm,~simstandard::get_model_implied_correlations(m=.,observed=TRUE,
                                                                            latent=FALSE,errors=FALSE))

  #simulate ordinal data and add d4 to scale back to original metric
  set.seed(269854)
  all_data_misspec <- purrr::map(Miss_Cor,~as.data.frame(GenOrd::ordsample(n=n*r, marginal=p1,Sigma=.)+d4))

  #calculate target matrix of fitted model
  #z<-GenOrd::ordcont(marginal=p1,Sigma=a)

  #get estimates from target matrix
  #norm_model<-lavaan::cfa(model=mod, sample.cov=z$SigmaC, sample.nobs=n)

  #Model statement for target matrix estimates
  #n_model<-cfa_lavmod(norm_model)

  #find misspecifications on target matrix scale
  #misspec_dgm <- DGM_Multi_HB(n_model)

  #get model-implied correlation for misspecification on target matr
  #all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*r,
    #                                                                        latent=FALSE,errors=FALSE))

  #n <- min(n,2000)

  #Set seed
  #set.seed(269854)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
 # r <- reps

  #all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*r,
      #                                                                      latent=FALSE,errors=FALSE))

  #Grab data level of the list
  #data <- purrr::map(misspec_data,2)

  #loop through misspecifications
  #for (x in 1:length(misspec_dgm))
 # {
  #  for (i in 1:ncol(dp))
   # {# first loop transforms highest category
    #  u<- as.numeric(sum(!is.na(dp[,i])))
     # all_data_misspec[[x]]<-all_data_misspec[[x]] %>%
      #  dplyr::mutate(!!a2[i] := case_when(
       #   !!rlang::sym(a2[i])  <= dp[1,i] ~100,
        #  !!rlang::sym(a2[i]) >   dp[u,i] ~100*(u+1),
         # TRUE ~ !!rlang::sym(a2[i]))
        #)
      #second loop transforms all lower categories
      #if(sum(!is.na(dp[,i])) > 1){
       # for (j in 1:sum(!is.na(dp[,i]))) {
        #  all_data_misspec[[x]]<-all_data_misspec[[x]] %>%
         #   dplyr::mutate(!!a2[i] := case_when(
          #    between( !!rlang::sym(a2[i]) , dp[j,i], dp[j+1,i]) ~ as.numeric(100*(j+1)),
           #   TRUE~ !!rlang::sym(a2[i]))
            #)
        #}
    #  }
  #  }
    #loops multiply by 100 to avoid overwriting MVN data with values above 1
    #divide by 100 to put things back onto Likert metric
   # all_data_misspec[[x]]<-all_data_misspec[[x]]/100
 # }

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- rep(1:r,n)

  #Combine indicator with dataset for each element in list
  dat_rep_misspec <- purrr::map(all_data_misspec,~cbind(.,rep_id_misspec))

  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~group_by(.,rep_id_misspec) %>%
                               tidyr::nest())

  #Grab data level of the list
  data_m <- purrr::map(misspec_data,2)

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa for each element in the list
  misspec_cfa <- purrr::map(data_m, function(x) purrr::map(x, function(y) lavaan::cfa(model=mod,
                                                                                      estimator=estimator,
                                                                                      data=y,
                                                                                      std.lv=TRUE,
                                                                                      se=se,
                                                                                      check.gradient=FALSE,
                                                                                      check.post=FALSE,
                                                                                      check.vcov=FALSE,
                                                                                      control=list(rel.tol=.001))))

  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, ind)) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))

  set.seed(NULL)

  return(misspec_fit_sum)

}

## Function for correct model

true_fit_multi_likert2 <- function(model,data,n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  true_dgm <- model

  factors <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(lhs != rhs)%>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #Get model-implied matrix of input model
  norm <- simstandard::sim_standardized_matrices(true_dgm)
  normx<-norm$Correlations$R
  l<-nrow(normx)-nrow(factors)
  a<-normx[1:l,1:l]
  diag(a)=1

  #names of variable in model
  names<-colnames(a)

  #identify names of variables from dataset
  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]

  #rescale so that minimum value is always ==1
  #needed to properly index computations below
  #also needed to avoid 0s because there is multiplication involved
  #will be scaled back at the end after loop indexes are not longer needed
  d2<-matrix(sapply(data1, function(x) min(x, na.rm=T)-1), nrow=1, ncol=ncol(data1))
  d3<-matrix(rep(d2,each=nrow(data1)), nrow=nrow(data1), ncol=ncol(data1))
  d4<-matrix(rep(d2,each=nrow(data1)), nrow=r*nrow(data1), ncol=ncol(data1))
  colnames(d3)<-colnames(data1)
  colnames(d4)<-colnames(data1)
  data1<-data1-d3

    #create empty matrix for proportions in each category (g)
  #GenOrd requires list (p) for each item
  #t is thresholds for eventual discretization
  g<-matrix(nrow=(max(data1,na.rm=T)-min(data1, na.rm=T)), ncol=ncol(data1))
  p<-list()
  t<-g

  #loop through all variables in fitted model
  for (i in 1:ncol(data1)){
    #loop through all categories for each variable
    for (h in min(data1[,i],na.rm=T):(max(data1[,i],na.rm=T)-1)){
      #proportion of responses at or below each category, only count non-missing in the denominator
      g[h,i]<-sum(table(data1[,i])[names(table(data1[,i]))<=h])/sum(!is.na(data1[,i]))
      p[[i]]<-g[,i]
      t[h,i]<-qnorm(c(g[h,i]))
      dp<-as.data.frame(t)
    }
  }

  #handle different number of categories per item
  p1<-list()
  for (x in 1:length(p)){
  p1[[x]]<-p[[x]][!is.na(p[[x]])]
  }

  #Assign column names so it's clear which thresholds go with which variable
  colnames(g)<-colnames(data1)
  colnames(dp)<-colnames(data1)
  a2<-colnames(g)


  #calculate target matrix of fitted model
  #z<-GenOrd::ordcont(marginal=p1,Sigma=a)

  #Set Seed
  set.seed(326267)
  #simultate ordinal data with same observed correlation matrix as observed data
  all_data_true<-as.data.frame(GenOrd::ordsample(n=n*r, marginal=p1, Sigma=a)+d4)
  #add d3 back to get original categories

  #get estimates from target matrix
  #norm_model<-lavaan::cfa(model=mod, sample.cov=z$SigmaC, sample.nobs=n)

  #Model statement for target matrix estimates
  #n_model<-cfa_lavmod(norm_model)

  #Use max sample size of 10000
  #n <- base::min(n,2000)

  #Set Seed
  #set.seed(326267)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  #r <- reps

  #Simulate one large dataset
  #all_data_true <- simstandard::sim_standardized(m=n_model,n = n*r,
  #                                               latent = FALSE,
  #                                               errors = FALSE)

 # for (i in 1:ncol(dp))
 #  {# first loop transforms highest category
  #  u<- as.numeric(sum(!is.na(dp[,i])))
   # all_data_true<-all_data_true %>%
    #  dplyr::mutate(!!a2[i] := case_when(
     #   !!rlang::sym(a2[i])  <= dp[1,i] ~100,
      #  !!rlang::sym(a2[i]) >   dp[u,i] ~100*(u+1),
       # TRUE ~ !!rlang::sym(a2[i]))
      #)
    #second loop transforms all lower categories
    #if(sum(!is.na(dp[,i])) > 1){
      #for (j in 1:sum(!is.na(dp[,i]))) {
        #all_data_true<-all_data_true %>%
          #dplyr::mutate(!!a2[i] := case_when(
            #between( !!rlang::sym(a2[i]) , dp[j,i], dp[j+1,i]) ~ as.numeric(100*(j+1)),
            #TRUE~ !!rlang::sym(a2[i]))
          #)
      #}
    #}
  #}
  #loops multiply by 100 to avoid overwriting MVN data with values above 1
  #divide by 100 to put things back onto Likert metric
  #all_data_true<-all_data_true/100

  rep_id_true <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)

  #Group and list
  data_t <- dat_rep_true %>%
    dplyr::group_by(rep_id_true) %>%
    tidyr::nest() %>%
    base::as.list()

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa
  true_cfa <- purrr::map(data_t$data,function(x) lavaan::cfa(model=mod,
                                                             estimator=estimator,
                                                             data=x,
                                                             std.lv=TRUE,
                                                             se=se,
                                                             check.gradient=FALSE,
                                                             check.post=FALSE,
                                                             check.vcov=FALSE,
                                                             control=list(rel.tol=.001)))

  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., ind)) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)
}

multi_df_likert2 <- function(model,data,n,estimator,reps){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- multi_fit_likert2(model,data,n,estimator,reps)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_multi_likert2(model,data,n,estimator,reps)

  #Produce final table by level
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}

data_likert2 <- function(model,data,n){

  mod <- cleanmodel(model)

  true_dgm <- model

  factors <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(lhs != rhs)%>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #Get model-implied matrix of input model
  norm <- simstandard::sim_standardized_matrices(true_dgm)
  normx<-norm$Correlations$R
  l<-nrow(normx)-nrow(factors)
  a<-normx[1:l,1:l]
  diag(a)=1

  #names of variable in model
  names<-colnames(a)

  #identify names of variables from dataset
  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]

  #rescale so that minimum value is always ==1
  #needed to properly index computations below
  #also needed to avoid 0s because there is multiplication involved
  #will be scaled back at the end after loop indexes are not longer needed
  d2<-matrix(sapply(data1, function(x) min(x, na.rm=T)-1), nrow=1, ncol=ncol(data1))
  d3<-matrix(rep(d2,each=nrow(data1)), nrow=nrow(data1), ncol=ncol(data1))
  #d4<-matrix(rep(d2,each=nrow(data1)), nrow=r*nrow(data1), ncol=ncol(data1))
  colnames(d3)<-colnames(data1)
  #colnames(d4)<-colnames(data1)
  data1<-data1-d3

  #create empty matrix for proportions in each category (g)
  #GenOrd requires list (p) for each item
  #t is thresholds for eventual discretization
  g<-matrix(nrow=(max(data1,na.rm=T)-min(data1, na.rm=T)), ncol=ncol(data1))
  p<-list()
  t<-g

  #loop through all variables in fitted model
  for (i in 1:ncol(data1)){
    #loop through all categories for each variable
    for (h in min(data1[,i],na.rm=T):(max(data1[,i],na.rm=T)-1)){
      #proportion of responses at or below each category, only count non-missing in the denominator
      g[h,i]<-sum(table(data1[,i])[names(table(data1[,i]))<=h])/sum(!is.na(data1[,i]))
      p[[i]]<-g[,i]
      t[h,i]<-qnorm(c(g[h,i]))
      dp<-as.data.frame(t)
    }
  }

  #handle different number of categories per item
  p1<-list()
  for (x in 1:length(p)){
    p1[[x]]<-p[[x]][!is.na(p[[x]])]
  }

  #Assign column names so it's clear which thresholds go with which variable
  colnames(g)<-colnames(data1)
  colnames(dp)<-colnames(data1)
  a2<-colnames(g)

  #calculate target matrix of fitted model
  #z<-GenOrd::ordcont(marginal=p1,Sigma=a)

  #Set Seed
  set.seed(326267)
  #simultate ordinal data with same observed correlation matrix as observed data
  all_data_true<-as.data.frame(GenOrd::ordsample(n=n, marginal=p1, Sigma=a)+d3)
  #add d3 back to get original categories

  miss_dat<-list()

  miss_dat$sim<-all_data_true
  miss_dat$orig<-data[,names]

  return(miss_dat)
}

##################################################
################# likertOne ######################
##################################################


### One-factor: Simulate fit indices for misspecified model for all levels ###

one_fit_likert <- function(model,data,n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm
  misspec_dgm <- DGM_one(model)

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set seed
  set.seed(649364)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Simulate one large dataset for each misspecification (use map to apply across each
  #element (set of misspecifications) in the list)
  all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*r,
                                                                            latent=FALSE,errors=FALSE))

  names<-colnames(all_data_misspec[[1]])

  #identify names of variables from dataset
  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]


  #create empty matrix for proportions (g) and threshold (p)
  g<-matrix(nrow=(max(data1)-min(data1)), ncol=ncol(data1))
  p<-g

  #loop through all variables in fitted model
  for (i in 1:ncol(data1)){
    #loop through all categories in fitted model
    for (h in min(data1[,i]):(max(data1[,i])-1)){
      #proportion of responses at or below each category
      g[h,i]<-sum(table(data1[,i])[names(table(data1[,i]))<=h])/nrow(data1)
      #inverse standard normal to determine thresholds
      p[h,i]<-qnorm(c(g[h,i]))
      dp<-as.data.frame(p)
    }
  }

  colnames(dp)<-colnames(data1)
  a2<-colnames(dp)


  #loop through misspecifications
  for (x in 1:length(misspec_dgm))
  {
    for (i in 1:ncol(dp))
    {# first loop transforms highest category
      u<- as.numeric(sum(!is.na(dp[,i])))
      all_data_misspec[[x]]<-all_data_misspec[[x]] %>%
        dplyr::mutate(!!a2[i] := case_when(
          !!rlang::sym(a2[i])  <= dp[1,i] ~100,
          !!rlang::sym(a2[i]) >   dp[u,i] ~100*(u+1),
          TRUE ~ !!rlang::sym(a2[i]))
        )
      #second loop transforms all lower categories
      if(sum(!is.na(dp[,i])) > 1){
        for (j in 1:sum(!is.na(dp[,i]))) {
          all_data_misspec[[x]]<-all_data_misspec[[x]] %>%
            dplyr::mutate(!!a2[i] := case_when(
              between( !!rlang::sym(a2[i]) , dp[j,i], dp[j+1,i]) ~ as.numeric(100*(j+1)),
              TRUE~ !!rlang::sym(a2[i]))
            )
        }
      }
    }
    #loops multiply by 100 to avoid overwriting MVN data with values above 1
    #divide by 100 to put things back onto Likert metric
    all_data_misspec[[x]]<-all_data_misspec[[x]]/100
  }


  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_misspec <- purrr::map(all_data_misspec,~base::cbind(.,rep_id_misspec))

  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~dplyr::group_by(.,rep_id_misspec) %>%
                               tidyr::nest())

  #Grab data level of the list
  data <- purrr::map(misspec_data,2)

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa for each element in the list
  misspec_cfa <- purrr::map(data, function(x) purrr::map(x, function(y) lavaan::cfa(model=mod,
                                                                                    estimator=estimator,
                                                                                    data=y,
                                                                                    std.lv=TRUE,
                                                                                    se=se,
                                                                                    check.gradient=FALSE,
                                                                                    check.post=FALSE,
                                                                                    check.vcov=FALSE,
                                                                                    control=list(rel.tol=.001))))

  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, ind)) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))

  set.seed(NULL)

  return(misspec_fit_sum)

}

#### One_Factor: Function to create True DGM (aka, just the model the user read in) ####

true_fit_one_likert <- function(model,data,n,estimator,reps){

  #Get clean model equation
  mod <- cleanmodel(model)

  true_dgm <- model

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set Seed
  set.seed(326267)

  #Number of reps (default is 500 and shouldn't be changed by empirical researchers)
  r <- reps

  #Simulate one large dataset
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n*r,
                                                 latent = FALSE,
                                                 errors = FALSE)

  names<-colnames(all_data_true)

  #identify names of variables from dataset
  data1<-data[,names]
  #remove cases if all relevant variables are NA
  data1<-data1[rowSums(is.na(data1)) != ncol(data1), ]

  #create empty matrix for proportions (g) and threshold (p)
  g<-matrix(nrow=(max(data1, na.rm=T)-min(data1, na.rm=T)), ncol=ncol(data1))
  p<-g

  #loop through all variables in fitted model
  for (i in 1:ncol(data1)){
    #loop through all categories in fitted model
    for (h in min(data1[,i], na.rm=T):(max(data1[,i], na.rm=T)-1)){
      #proportion of responses at or below each category
      g[h,i]<-sum(table(data1[,i])[names(table(data1[,i]))<=h])/nrow(data1)
      #inverse standard normal to determine thresholds
      p[h,i]<-qnorm(c(g[h,i]))
      dp<-as.data.frame(p)
    }
  }

  colnames(dp)<-colnames(data1)
  a2<-colnames(dp)

  for (i in 1:ncol(dp))
  {# first loop transforms highest category
    u<- as.numeric(sum(!is.na(dp[,i])))
    all_data_true<-all_data_true %>%
      dplyr::mutate(!!a2[i] := case_when(
        !!rlang::sym(a2[i])  <= dp[1,i] ~100,
        !!rlang::sym(a2[i]) >   dp[u,i] ~100*(u+1),
        TRUE ~ !!rlang::sym(a2[i]))
      )
    #second loop transforms all lower categories
    if(sum(!is.na(dp[,i])) > 1){
      for (j in 1:sum(!is.na(dp[,i]))) {
        all_data_true<-all_data_true %>%
          dplyr::mutate(!!a2[i] := case_when(
            between( !!rlang::sym(a2[i]) , dp[j,i], dp[j+1,i]) ~ as.numeric(100*(j+1)),
            TRUE~ !!rlang::sym(a2[i]))
          )
      }
    }
  }

  #loops multiply by 100 to avoid overwriting MVN data with values above 1
  #divide by 100 to put things back onto Likert metric
  all_data_true<-all_data_true/100

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_true <- base::rep(1:r,n)

  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)

  #Group and list
  true_data <- dat_rep_true %>%
    dplyr::group_by(rep_id_true) %>%
    tidyr::nest() %>%
    base::as.list()

  #if using ULS or one of its variants, the following must be specified for standard errors: se="standard"
  se<-ifelse(startsWith(estimator,"ULS"),"standard","none")
  if(endsWith(estimator,"MVS") | endsWith(estimator,"V")) {
    ind<-c("srmr","rmsea.scaled","cfi.scaled")
  } else if(endsWith(estimator,"M") | endsWith(estimator,"R")) {
    ind<-c("srmr","rmsea.robust","cfi.robust")
  }  else {
    ind<-c("srmr","rmsea","cfi")
  }

  #Run 500 cfa
  true_cfa <- purrr::map(true_data$data,function(x) lavaan::cfa(model=mod,
                                                                estimator=estimator,
                                                                data=x,
                                                                std.lv=TRUE,
                                                                se=se,
                                                                check.gradient=FALSE,
                                                                check.post=FALSE,
                                                                check.vcov=FALSE,
                                                                control=list(rel.tol=.001)))

  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., ind)) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#### One-Factor: Function to combine both model fit stats for all levels into one dataframe ####

one_df_likert <- function(model,data,n,estimator,reps){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- one_fit_likert(model,data,n,estimator,reps)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_one_likert(model,data,n,estimator,reps)

  #Produce final table by level
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}

##################################################
################# DDDFI ##########################
##################################################

#Function to strip estimates from model statement
#Used to apply fitted model to simulated data
#different from original cleanmodel function to allow..
#... for broader types of paths
cleanmodel_3DFI <- function(model){

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


#There is no build in function for the LKJ distribution
#Manually write out LKJ distribution (taken trialr package, maintained by Kristian Brock)
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

#Creates discrepancy matrices
#generates misspecified data
#fits original model to misspecified data
#collates fit indices into a data frame
miss_fit <- function(model,data,n,reps,estimator,MAD,scale){

  #strip estimates from model statement
  mod <- cleanmodel_3DFI(model)

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
  n <- n

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

#generate data consistent with the original model
#fit orginal model to correct/consistent data
#collated fit indices into a data frame

true_fit<- function(model,data,n,reps, estimator, MAD,scale){

  #strip estimates from model statement
  mod <- cleanmodel_3DFI(model)

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
  n <-n

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

#Function to combine all indices into one dataframe
#Will used to created distributions that determine optimal cutoff
combined <- function(model,data,n,reps,estimator,MAD,scale){

  #Use max sample size of 5000
  n <- n

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
