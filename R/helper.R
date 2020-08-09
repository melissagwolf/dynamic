#' @importFrom dplyr filter select mutate full_join group_by add_tally
#' ungroup slice arrange summarise count as_tibble pull recode
#' @importFrom tidyr pivot_longer unite
#' @importFrom lavaan lavaanify cfa fitMeasures
#' @importFrom simstandard sim_standardized
#' @import ggplot2
#' @import magrittr
#' @importFrom stats quantile
#' @noRd

#### Function for Number of Factors ####

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
unstandardized <- function(model){

  lav_file <- lavaan::lavaanify(model, fixed.x=FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs)

  one_plus <- lav_file %>%
    dplyr::filter(ustart >= 1) %>%
    base::nrow()

  return(one_plus)
}



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

#### Function to create model statement without numbers from user model (for input) ####

cleanmodel <- function(model){

  model %>%
    lavaan::lavaanify(fixed.x = FALSE) %>%
    dplyr::filter(.data$lhs != .data$rhs) %>%
    dplyr::group_by(.data$lhs, .data$op) %>%
    dplyr::summarise(rhs = paste(.data$rhs, collapse = " + ")) %>%
    dplyr::arrange(dplyr::desc(.data$op)) %>%
    tidyr::unite("l", .data$lhs, .data$op, .data$rhs, sep = " ") %>%
    dplyr::pull(.data$l)

}

#Catch regular warning "some estimated ov variances are negative"
#Use in misspecified_model_fit function with cfa
#http://romainfrancois.blog.free.fr/index.php?post/2009/05/20/Disable-specific-warnings

hide_ov <- function(h){
  if(any(grepl("some estimated ov variances are negative", h)))
    invokeRestart("muffleWarning")
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
