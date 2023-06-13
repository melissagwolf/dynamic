devtools::install_github("melissagwolf/Dynamic",  force=T)

library(dplyr)
dynamic::

  head(data)
data<-read.csv(file.choose())

m<-"E=~ bfi_e1+bfi_e2+bfi_e3+bfi_e4+bfi_e5+bfi_e6+bfi_e7+bfi_e8"

model<-lavaan::cfa(m, data=data)
lavaan::summary(m2)
n<-500
reps<-5
estimator<-"ML"
model2<-cfa_lavmod(m2)


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
      g[h,i]<-sum(tabulate(data1[,i])[names(table(data1[,i]))<=h])/nrow(data1)
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
  g<-matrix(nrow=(max(data1)-min(data1)), ncol=ncol(data1))
  p<-g

  #loop through all variables in fitted model
  for (i in 1:ncol(data1)){
    #loop through all categories in fitted model
    for (h in min(data1[,i]):(max(data1[,i])-1)){
      #proportion of responses at or below each category
      g[h,i]<-sum(tabulate(data1[,i])[names(table(data1[,i]))<=h])/nrow(data1)
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

likertOne(m1, data=data, rep=5)

cleanmodel(model2)

devtools::install_github("melissagwolf/dynamic", force=T)

library(devtools)
devtools::document()


dynamic::l
