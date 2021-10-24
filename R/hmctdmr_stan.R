hmctdmr_stan_container <- function(object){
  stan_model <- object %>% hmctdmr_stan_model()

  stanD <- hmctdmr_stan_stanD(object$data, object$prior)
  stanInit <- hmctdmr_stan_stanInit(object$drug, object$stan_sample_option)

  fit <- hmctdmr_stan_sample(object, stan_model, stanD, stanInit)

  stanOut <- list(
    stan_model=stan_model,
    stan_train=stanD,
    stan_init=stanInit,
    fit=fit
  )
  return(stanOut)
}

hmctdmr_stan_model <- function(object){

  if(is.null(object$stan_model_option$stan_file)){
    object$stan_model_option$stan_file <- search_stan_file(object$drug)
  }

  stan_mod <- do.call(cmdstanr::cmdstan_model, object$stan_model_option)

  return(stan_mod)
}

hmctdmr_stan_stanD <- function(data, prior){

  nt <- nrow(data)
  iObs <- with(data, (1:nrow(data))[evid == 0])
  nObs <- length(iObs)
  cObs <- with(data, DV[iObs])
  return (c((with(data, list( 
                      nt=nt, nObs=nObs,
                      iObs=as.array(iObs), cObs=as.array(cObs), 
                      cmt=as.array(cmt), time=as.array(time), 
                      amt=as.array(amt), rate= as.array(rate), 
                      evid=as.array(evid), ii=as.array(ii), 
                      addl=as.array(addl), ss=as.array(ss),
                      WT=as.array(WT), LBW=as.array(LBW),
                      CrCL=as.array(CrCL)))), prior))  
}

hmctdmr_stan_stanInit <- function(drug, stan_sample_option=list(...)){

  nChains <- stan_sample_option$chains
  create.stanInit <- switch(drug,
                "amikacin" = function(chains){
                              list( CL_NR=0.0417,
                                    CL_SLOPE=0.815,
                                    VD_NR = 0.27)
                              },
                "vancomycin" = function(chains){
                              list( CL_NR=0.0417,
                                    CL_SLOPE=0.815,
                                    VD_NR = 0.27)
                              },
                "theophiline" = function(chains){
                              list( CL_NR=0.0417,
                                    CL_SLOPE=0.815,
                                    VD_NR = 0.27)
                              },
                "phenytoin" = function(chains){
                              list(ln_Vmax=log(500),  
                                    ln_Km=log(5.0),
                                    ln_VD_NR=log(0.8))
                              },
                            NULL
                          )

  stanInit <- lapply(1:nChains, create.stanInit)

  return(stanInit)
}

hmctdmr_stan_sample <- function(object, stan_model, stanD, stanInit=NULL){

  object$stan_sample_option$data <- stanD

  if( !("refresh" %in% names(object$stan_sample_option)) )
    object$stan_sample_option$refresh <- 0

  if(!is.null(stanInit)) 
    object$stan_sample_option$init <- stanInit
  
  fit <- suppressMessages(do.call(stan_model$sample, object$stan_sample_option))
  # fit <- (do.call(stan_model$sample, object$stan_sample_option))

  return(fit)

}

search_stan_file <- function(drug){
  message("search stan_file ..")
  stan_file <- switch(drug,
                "amikacin" = file.path(stanLib(), "amikacin.stan"),
                "vancomycin" = file.path(stanLib(), "vancomycin.stan"),
                "theophiline" = file.path(stanLib(), "theophiline.stan"),
                "phenytoin" = file.path(stanLib(), "phenytoin.stan"),
                NULL
              )
  if(!is.null(stan_file)) { message("find stan_file : ", stan_file) }
  return(stan_file)
}

match_stanD <- function(drug){
  stanD <- switch(drug,
                "amikacin" = function(){message("STAND::AMI")},
                "vancomycin" = file.path(stanLib(), "vancomycin.stan"),
                "theophiline" = file.path(stanLib(), "theophiline.stan"),
                "phenytoin" = file.path(stanLib(), "phenytoin.stan"),
                NULL
              )

}