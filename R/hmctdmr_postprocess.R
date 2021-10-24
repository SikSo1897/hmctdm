hmctdmr_postprocess <- function(hmctdm){

  post <- list()
  post$stan_model_option <- cbind(
                              option="stan_model_option",
                              do.call(cbind, hmctdm$stan_model_option), row.names=NULL
                            ) %>% as_tibble

  
  post$stan_sample_option <- cbind(
                              option="stan_sample_option",
                              do.call(cbind, hmctdm$stan_sample_option), row.names=NULL
                            ) %>% as_tibble

  fit <- hmctdm$stan_result$fit$summary()

  not_pk_pattern <- "cHat|lp|theta|x|sigma"

  pk_parm <- fit %>% filter(!stringr::str_detect(variable, not_pk_pattern))
  
  idx.cHat <- (1:nrow(fit))[stringr::str_detect(fit$variable, "cHat\\[")]

  cHat <- fit[idx.cHat, "median"] %>% rename(cHat=median)
  
  post$status <- hmctdm$status
  post$drug <- hmctdm$drug
  post$stan_model_option
  post$stan_sample_option
  post$param <- pk_parm
  post$data <- hmctdm$data
  post$est <- hmctdm$data
  post$est$cHat <- with(cHat, as.array(cHat))

  print(post$data)

  idata <- post$data %>% 
            filter(evid==1)  %>%
            select(any_of(c("ID","SEX","AGE","WT","LBW","HT","SCR","CrCL"))) %>%
            cbind(
              do.call(cbind, setNames(as.list(pk_parm$median), pk_parm$variable))
            )
  
  post$idata <- idata

  class(post) <- "hmctdm"

  return(post)
  
}