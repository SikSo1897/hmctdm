hmctdmr_stan_container <- function(object){
  stan_model <- object %>% hmctdmr_stan_model()
  # print(object$data)

  # print(object$stan_sample_option)
  stanD <- hmctdmr_stan_stanD(object$drug, object$data)
  stanInit <- hmctdmr_stan_stanInit(object$drug, object$stan_sample_option)
  # print(stanInit)

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

  # print(object$stan_model_option)
  # print(unlist(object$stan_model_option, use.names=FALSE))
  # print(object$stan_model_option)

  stan_mod <- do.call(cmdstanr::cmdstan_model, object$stan_model_option)

  return(stan_mod)
}

hmctdmr_stan_stanD <- function(drug, data){
  get.stanD <- switch(drug,
    "amikacin" = function(data=data){ 
                  nt <- nrow(data)
                  iObs <- with(data, (1:nrow(data))[evid == 0])
                  nObs <- length(iObs)
                  cObs <- with(data, DV[iObs])
                  # time <- with(data, time[time != 0])

                  return (with(data, list( 
                                      nt=nt, nObs=nObs,
                                      iObs=as.array(iObs), cObs=as.array(cObs), 
                                      cmt=as.array(cmt), time=as.array(time), 
                                      amt=as.array(amt), rate= as.array(rate), 
                                      evid=as.array(evid), ii=as.array(ii), 
                                      addl=as.array(addl), ss=as.array(ss), 
                                      CrCL=as.array(CrCL), LBW=as.array(LBW))) )
                },
    "vancomycin" = file.path(stanLib(), "vancomycin.stan"),
    "theophiline" = file.path(stanLib(), "theophiline.stan"),
    "phenytoin" = file.path(stanLib(), "phenytoin.stan"),
                NULL
              )
  stanD <- get.stanD(data)

  return(stanD)
  
}

hmctdmr_stan_stanInit <- function(drug, stan_sample_option=list(...)){
  nChains <- stan_sample_option$chains
  get.stanInit <- switch(drug,
                "amikacin" = function(chains){
                              list( CL_NR=0.0417,
                                    CL_SLOPE=0.815,
                                    VD_NR = 0.27)
                              },
                "vancomycin" = file.path(stanLib(), "vancomycin.stan"),
                "theophiline" = file.path(stanLib(), "theophiline.stan"),
                "phenytoin" = file.path(stanLib(), "phenytoin.stan"),
                            NULL
                          )

  stanInit <- lapply(1:nChains, get.stanInit)

  return(stanInit)
}

hmctdmr_stan_sample <- function(object, stan_model, stanD, stanInit=NULL){

  object$stan_sample_option$data <- stanD

  if( !("refresh" %in% names(object$stan_sample_option)) )
    object$stan_sample_option$refresh <- 0

  if(!is.null(stanInit)) 
    object$stan_sample_option$init <- stanInit
  
  fit <- suppressMessages(do.call(stan_model$sample, object$stan_sample_option))
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