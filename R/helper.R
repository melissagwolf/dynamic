#' @importFrom dplyr filter select mutate full_join group_by add_tally
#' ungroup slice arrange summarise count as_tibble pull recode distinct
#' mutate_if slice_max slice_min
#' @importFrom tidyr pivot_longer unite nest extract
#' @importFrom lavaan lavaanify cfa fitMeasures
#' @importFrom simstandard sim_standardized
#' @importFrom purrr map map_dfr
#' @import ggplot2
#' @importFrom stats quantile


############################################################################
########################### Multiple Functions #############################
############################################################################

#### Function to create model statement without numbers from user model (for input) ####
#Copy from OG
cleanmodel <- function(model){

  suppressMessages(model %>%
                     lavaan::lavaanify(fixed.x = FALSE) %>%
                     dplyr::filter(.data$lhs != .data$rhs) %>%
                     dplyr::group_by(.data$lhs, .data$op) %>%
                     dplyr::summarise(rhs = paste(.data$rhs, collapse = " + ")) %>%
                     dplyr::arrange(dplyr::desc(.data$op)) %>%
                     tidyr::unite("l", .data$lhs, .data$op, .data$rhs, sep = " ") %>%
                     dplyr::pull(.data$l))

}

#### Function for Number of Factors ####
#Copy from OG

number_factor <- function(model){

  #prep the model
  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)

  #isolate factors
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #identify number of factors in model
  num_factors <- base::nrow(factors)

  return(num_factors)
}

#Did they enter unstandardized loadings?  Aka, do they have any loadings = 1?
#Copy from OG
#Used for error message

unstandardized <- function(model){

  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)

  one_plus <- lav_file %>%
    dplyr::filter(ustart >= 1) %>%
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

#Catch regular warning "some estimated ov variances are negative"
#Use in misspecified_model_fit function with cfa
#http://romainfrancois.blog.free.fr/index.php?post/2009/05/20/Disable-specific-warnings

hide_ov <- function(h){
  if(any(grepl("some estimated ov variances are negative", h)))
    invokeRestart("muffleWarning")
}

hide_objtest <- function(j){
  if(any(grepl("restarting interrupted promise evaluation", j)))
    invokeRestart("muffleWarning")
}

##### NEW: Extract n from lavaan object #####
#SPECIFIC TO R PACKAGE

cfa_n <- function(model){

  #Extract n from lavaan object
  #Warning message to hide warning when someone enters a non-lavaan object (error message will display instead)
  #Only need it here because this is the first argument in cfaHB and cfaOne
  n <- base::withCallingHandlers(base::unlist(model@SampleStats@nobs),warning=hide_objtest)
  return(n)
}

##### NEW: Extract model statement from lavaan object #####
#SPECIFIC TO R PACKAGE

cfa_lavmod <- function(model){

  #Extract standardized solution from lavaan object
  lav <- lavaan::standardizedSolution(model)

  #Create model statement
  ss_mod <- suppressMessages(lav %>%
                               dplyr::filter(lhs != rhs) %>%
                               dplyr::group_by(lhs,op) %>%
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

############################################################################
################################# cfaFit ###################################
############################################################################


#### Function for Single Factor Misspecification (Correlation) ####

single_factor <- function(model){

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)

  #identify the factor name
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
    dplyr::ungroup()

  #Select the item with the lowest loading
  lowest1 <- solo_items %>%
    dplyr::arrange(abs(ustart)) %>%
    dplyr::slice(1)

  #Select the item with the second lowest loading
  lowest2 <- solo_items %>%
    dplyr::arrange(abs(ustart)) %>%
    dplyr::slice(2) %>%
    `colnames<-`(c("lhs2","op2","rhs2","ustart2","n2"))

  #Compute the model implied residual correlation (simstandard)
  U1 <- lowest1$ustart
  U2 <- lowest2$ustart2
  Cor <- U1*U2
  #U1E <- (1-(U1^2))              #for dg by lavaan
  #U2E <- (1-(U2^2))
  #Cov <- Cor*sqrt(U1E*U2E)

  #Create model DF
  Residual_Correlation <- Cor %>%
    base::as.data.frame() %>%
    `colnames<-`("Cor") %>%
    dplyr::mutate(Cor=round(Cor,4)) %>%
    base::cbind(lowest1,lowest2) %>%
    dplyr::mutate(operator="~~") %>%
    dplyr::mutate(times="*") %>%
    dplyr::select(rhs,operator,Cor,times,rhs2) %>%
    tidyr::unite("V1",sep=" ")

  return(Residual_Correlation)
}

#### Function for multi-factor misspecification (Cross-loading) ####

multi_factor <- function(model){

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)

  #identify all factor names
  factors <- lav_file %>%
    dplyr::filter(op=="=~") %>%
    dplyr::select(lhs) %>%
    base::unique()

  #Compute Coefficient H for each factor
  Coef_H <- lavaan::lavaanify(model, fixed.x = FALSE) %>%
    dplyr::filter(lhs != rhs) %>%
    dplyr::filter(op == "=~") %>%
    dplyr::mutate(L_Sq=ustart^2) %>%
    dplyr::mutate(E_Var=1-L_Sq) %>%
    dplyr::mutate(Div=L_Sq/E_Var) %>%
    dplyr::group_by(lhs) %>%
    dplyr::summarise(Sum=sum(Div)) %>%
    dplyr::mutate(H=((1+(Sum^-1))^-1)) %>%
    dplyr::select(-Sum) %>%
    dplyr::arrange(-H) %>%
    `colnames<-`(c("rhs","H"))

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
    base::as.data.frame() %>%
    `colnames<-`(c("lhs","Remaining"))

  #Figure out which items are ideally eligible for cross-loading misspecification
  eligible <- remaining %>%
    dplyr::full_join(num_items,by="lhs") %>%
    dplyr::filter(Original>2 & Remaining != "NA") %>%
    dplyr::mutate(Eligible=1) %>%
    dplyr::select(lhs,Eligible)

  #Compute number that are ideally eligible
  Num_Eligible <- base::nrow(eligible)

  #Identify item and factor with lowest loading (that doesn't have an existing error cov)
  #We prefer to select a factor with 3+ items for identification purposes
  #If there is an eligible factor with at least 3 items, pick the lowest loading from one of those
  #If there isn't, then you can still proceed
  if(Num_Eligible>0){
    crosses <- solo_items %>%
      dplyr::filter(op=="=~") %>%
      dplyr::select(lhs,op,rhs,ustart) %>%
      dplyr::group_by(rhs) %>%                   #this is where we remove
      dplyr::add_tally() %>%                     #any items that load on
      dplyr::filter(n==1) %>%                    #more than one factor
      dplyr::full_join(eligible,by="lhs") %>%    #add factor eligibility status
      dplyr::filter(Eligible == 1) %>%           #select only factors that are eligible
      dplyr::arrange(abs(ustart)) %>%            #choose the lowest loading item
      base::as.data.frame() %>%                  #from a factor with more than 2 items
      dplyr::slice(1) %>%
      dplyr::select(-n,-op, -Eligible) %>%
      dplyr::as_tibble() %>%
      `colnames<-`(c("Factor","Item","Loading"))
  } else {
    crosses <- solo_items %>%
      dplyr::filter(op=="=~") %>%
      dplyr::select(lhs,op,rhs,ustart) %>%
      dplyr::group_by(rhs) %>%                   #this is where we remove
      dplyr::add_tally() %>%                     #any items that load on
      dplyr::filter(n==1) %>%                    #more than one factor
      dplyr::arrange(abs(ustart)) %>%            #choose the lowest loading item
      base::as.data.frame() %>%                 #from any factor
      dplyr::slice(1) %>%
      dplyr::select(-n,-op) %>%
      dplyr::as_tibble() %>%
      `colnames<-`(c("Factor","Item","Loading"))
  }

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

  #combine and clean
  modinfo <- factcor1 %>%
    dplyr::full_join(factcor2, by = c("lhs", "op", "rhs", "ustart", "type")) %>%
    base::cbind(crosses) %>%
    dplyr::full_join(Coef_H,by="rhs") %>%
    dplyr::filter(lhs == Factor) %>%
    dplyr::arrange(-H) %>%
    dplyr::slice(1) %>%
    dplyr::mutate(operator="=~")

  #Compute maximum allowable cross loading value
  F1 <- modinfo$ustart
  F1_Sq <- F1^2                       #.95 = monte carlo correction
  L1 <- modinfo$Loading               #to ensure that the dgm matrix
  L1_Sq <- L1^2                       #is positive definite
  E <- 1-L1_Sq
  MaxAllow <- ((sqrt((L1_Sq*F1_Sq)+E)-(L1*F1))*.95)

  #extract value of loading
  Final_Loading <- round(min(L1,MaxAllow),4)

  #Create model DF
  Cross_Loading <- modinfo %>%
    dplyr::select(rhs,Item,operator) %>%
    base::cbind(Final_Loading) %>%
    dplyr::mutate(times="*") %>%
    dplyr::select(rhs,operator,Final_Loading,times,Item) %>%
    tidyr::unite("V1",sep=" ")

  #return value to append to model statement
  return(Cross_Loading)
}

#### Function to create Misspecified DGM given the number of factors ####

Misspecified_DGM <- function(model){

  factor <- number_factor(model)

  if (factor > 1){
    multi <- multi_factor(model)
    multi_mod <- base::rbind(model,multi)
    multi_mod_c <- multi_mod$V1
    return(multi_mod_c)

  } else{
    single <- single_factor(model)
    single_mod <- base::rbind(model,single)
    single_mod_c <- single_mod$V1
    return(single_mod_c)
  }

}

#Simulate misspecified data fit stats

misspecified_model_fit <- function(model,n){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm
  misspec_dgm <- Misspecified_DGM(model)

  #Create empty df to put fit stats in
  misspec_fit <- data.frame(matrix(nrow=500,ncol=3))

  #Use max sample size of 10000
  n <- min(n,10000)

  #Simulate data through loop 500 times
  set.seed(649364)
  for (i in 1:500){
    misspec_data <- simstandard::sim_standardized(m=misspec_dgm,n = n,
                                                  latent = FALSE,
                                                  errors = FALSE)
    misspec_cfa <- base::withCallingHandlers(
      lavaan::cfa(model = mod, data = misspec_data, std.lv=TRUE), warning = hide_ov)
    misspec_fits <- lavaan::fitMeasures(misspec_cfa, c("srmr","rmsea","cfi"))
    misspec_fit[i,] <- misspec_fits
  }
  set.seed(NULL)

  #Clean up data
  misspec_fit_sum <- misspec_fit %>%
    `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
    dplyr::mutate(Type_M="Misspecified")

  return(misspec_fit_sum)

}

#Simulate true data fit stats

true_model_fit <- function(model,n){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for true dgm
  true_dgm <- model

  #Create empty df to put fit stats in
  true_fit <- data.frame(matrix(nrow=500,ncol=3))

  #Use max sample size of 10000
  n <- min(n,10000)

  #Simulate data through loop 500 times
  set.seed(326267)
  for (i in 1:500){
    true_data <- simstandard::sim_standardized(m=true_dgm,n = n,
                                               latent = FALSE,
                                               errors = FALSE)
    true_cfa <- base::withCallingHandlers(
      lavaan::cfa(model = mod, data = true_data, std.lv=TRUE), warning=hide_ov)
    true_fits <- lavaan::fitMeasures(true_cfa, c("srmr","rmsea","cfi"))
    true_fit[i,] <- true_fits
  }
  set.seed(NULL)

  #Clean up data
  true_fit_sum <- true_fit %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  return(true_fit_sum)

}

#Combine into one dataframe
dynamic_fit <- function(model,n){

  #Use max sample size of 10000
  n <- min(n,10000)

  #Get fit stats for misspecified model
  misspec_fit <- misspecified_model_fit(model,n)

  #Get fit stats for correctly specified model
  true_fit <- true_model_fit(model,n)

  #ifelse statements to produce final table
  Table <- base::cbind(misspec_fit,true_fit)

  #Final table
  return(Table)
}



############################################################################
################################# cfaOne ###################################
############################################################################



### One-factor: Function to see which items are available ###
## This name is new!!!

one_num <- function(model){

  #Rename (just to be consistent with shiny app)
  Mod_C <- model

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(Mod_C, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)

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
#This function name is new!!!!

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
## This name is new!!!

one_fit <- function(model,n){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm
  misspec_dgm <- DGM_one(model)

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set seed
  set.seed(649364)

  #Simulate one large dataset for each misspecification (use map to apply across each
  #element (set of misspecifications) in the list)
  all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*500,
                                                                            latent=FALSE,errors=FALSE))

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- base::rep(1:500,n)

  #Combine indicator with dataset
  dat_rep_misspec <- purrr::map(all_data_misspec,~base::cbind(.,rep_id_misspec))

  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~dplyr::group_by(.,rep_id_misspec) %>%
                               tidyr::nest())

  #Grab data level of the list
  data <- purrr::map(misspec_data,2)

  #Run 500 cfa
  misspec_cfa <- purrr::map(data, function(x) purrr::map(x, function(y) base::withCallingHandlers(lavaan::cfa(model = mod, data=y, std.lv=TRUE),
                                                                                                  warning=hide_ov)))

  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, c("srmr","rmsea","cfi"))) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))

  set.seed(NULL)

  return(misspec_fit_sum)

}

#### One_Factor: Function to create True DGM (aka, just the model the user read in) ####
## This name is new!!

true_fit_one <- function(model,n){

  #Get clean model equation
  mod <- cleanmodel(model)

  true_dgm <- model

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set Seed
  set.seed(326267)

  #Simulate one large dataset
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n*500,
                                                 latent = FALSE,
                                                 errors = FALSE)

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_true <- base::rep(1:500,n)

  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)

  #Group and list
  true_data <- dat_rep_true %>%
    dplyr::group_by(rep_id_true) %>%
    tidyr::nest() %>%
    base::as.list()

  #Run 500 cfa
  true_cfa <- purrr::map(true_data$data,~base::withCallingHandlers(lavaan::cfa(model = mod, data=., std.lv=TRUE),
                                                                   warning=hide_ov))

  #Extract fit stats from each rep (list) into a data frame and clean
  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., c("srmr","rmsea","cfi"))) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#### One-Factor: Function to combine both model fit stats for all levels into one dataframe ####
## New name!!

one_df <- function(model,n){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- one_fit(model,n)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_one(model,n)

  #Produce final table by level
  Table <- purrr::map(misspec_fit,~cbind(.,true_fit))

  #Final table
  return(Table)
}


############################################################################
################################# cfaHB ####################################
############################################################################


### Multi-factor: Function to see which items are available ###
## This name is new!!!

multi_num_HB <- function(model){

  #Rename (just to be consistent with shiny app)
  Mod_C <- model

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(Mod_C, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)

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
#This function name is new!!!!

multi_add_HB <- function(model){

  #read in the model
  Mod_C <- model

  #Lavaanify it - have lavaan tell us the parameters
  lav_file <- lavaan::lavaanify(Mod_C, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)

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
                     dplyr::summarise(Sum=sum(Div)) %>%
                     dplyr::mutate(H=((1+(Sum^-1))^-1)) %>%
                     dplyr::select(-Sum) %>%
                     dplyr::arrange(-H) %>%
                     `colnames<-`(c("rhs","H")))

  #isolate factors and factor correlations
  factcor1 <- factors %>%
    dplyr::mutate(type="Factor") %>%
    dplyr::full_join(lav_file, by = "lhs") %>%
    dplyr::mutate(type=recode(type, .missing ="Error Correlation")) %>%
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
    dplyr::mutate(MaxAllow=((base::sqrt(((L1_Sq*F1_Sq)+E))-(L1*F1))*.95),
                  MaxAllow2=base::round(MaxAllow,digits=4),
                  Final_Loading=base::pmin(Loading,MaxAllow2),
                  times="*") %>%
    dplyr::select(rhs,operator,Final_Loading,times,Item) %>%
    tidyr::unite("V1",sep=" ")

  #return value to append to model statement
  return(Cross_Loading)
}


#### Multi-factor: Function to create Misspecified DGM given the number of factors ####
## This name is new!!!

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
## This name is new!!!

multi_fit_HB <- function(model,n){

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for misspecified dgm (this is a list)
  misspec_dgm <- DGM_Multi_HB(model)

  #Use max sample size of 2000
  n <- min(n,2000)

  #Set seed
  set.seed(269854)

  #Simulate one large dataset for each misspecification (use map to apply across each
  #element (set of misspecifications) in the list)
  all_data_misspec <- purrr::map(misspec_dgm,~simstandard::sim_standardized(m=.,n=n*500,
                                                                            latent=FALSE,errors=FALSE))

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_misspec <- rep(1:500,n)

  #Combine indicator with dataset for each element in list
  dat_rep_misspec <- purrr::map(all_data_misspec,~cbind(.,rep_id_misspec))

  #Group and list
  misspec_data <- purrr::map(dat_rep_misspec,~group_by(.,rep_id_misspec) %>%
                               tidyr::nest())

  #Grab data level of the list
  data <- purrr::map(misspec_data,2)

  #Run 500 cfa for each element in the list
  misspec_cfa <- purrr::map(data, function(x) purrr::map(x, function(y) base::withCallingHandlers(lavaan::cfa(model = mod, data=y, std.lv=TRUE),
                                                                                                  warning = hide_ov)))

  #Extract fit stats from each rep (list) into a data frame and clean using nested lapply
  #map_dfr returns data frame instead of list
  #for each misspecification level (in the list), access the lavaan objects (x)
  #and extract the fit stats (y) - and return as a df
  misspec_fit_sum <- purrr::map(misspec_cfa, function(x) purrr::map_dfr(x, function(y) lavaan::fitMeasures(y, c("srmr","rmsea","cfi"))) %>%
                                  `colnames<-`(c("SRMR_M","RMSEA_M","CFI_M")) %>%
                                  dplyr::mutate(Type_M="Misspecified"))

  set.seed(NULL)

  return(misspec_fit_sum)

}

#### Multi_Factor: Function to create True DGM (aka, just the model the user read in) ####
## This name is new!!

true_fit_HB <- function(model,n){

  #Can make this faster by only doing it once
  #Would need to change table. Not sure what would happen to plot.
  #Already did this

  #Get clean model equation
  mod <- cleanmodel(model)

  #Get parameters for true DGM
  true_dgm <- model

  #Use max sample size of 10000
  n <- base::min(n,2000)

  #Set Seed
  set.seed(267326)

  #Simulate one large dataset
  all_data_true <- simstandard::sim_standardized(m=true_dgm,n = n*500,
                                                 latent = FALSE,
                                                 errors = FALSE)

  #Create indicator to split into 500 datasets for 500 reps
  rep_id_true <- base::rep(1:500,n)

  #Combine indicator with dataset
  dat_rep_true <- base::cbind(all_data_true,rep_id_true)

  #Group and list
  true_data <- dat_rep_true %>%
    dplyr::group_by(rep_id_true) %>%
    tidyr::nest() %>%
    base::as.list()

  #Run 500 cfa
  true_cfa <- purrr::map(true_data$data,~base::withCallingHandlers(lavaan::cfa(model = mod, data=., std.lv=TRUE),
                                                                   warning=hide_ov))

  #Extract fit stats from each rep (list) into a data frame and clean
  true_fit_sum <- purrr::map_dfr(true_cfa,~lavaan::fitMeasures(., c("srmr","rmsea","cfi"))) %>%
    `colnames<-`(c("SRMR_T","RMSEA_T","CFI_T")) %>%
    dplyr::mutate(Type_T="True")

  set.seed(NULL)

  return(true_fit_sum)

}

#### Multi-Factor: Function to combine both model fit stats for all levels into one dataframe ####
## New name!!

multi_df_HB <- function(model,n){

  #Use max sample size of 2000
  n <- min(n,2000)

  #Get fit stats for misspecified model
  misspec_fit <- multi_fit_HB(model,n)

  #Get fit stats for correctly specified model
  true_fit <- true_fit_HB(model,n)

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

  z=qnorm(1-alpha)
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
  T_ml <- base::withCallingHandlers(base::unlist(obj@test[["standard"]][["stat"]]), warning=hide_objtest)
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
