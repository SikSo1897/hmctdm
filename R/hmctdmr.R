#' hmctdmr() Function
#'
#' This function adds a and b
#' @param a ,b : two numbers to be operated
#' @keywords hmctdmr
#' @export
#' @examples
#' hmctdmr()
# _hmctdmr <- function(
#   data, 

# ){}

# HOME_TORSTAN <- "/home/sikso/work/00_Private/simulation/lib/Torsten"
# suppressMessages(set_cmdstan_path(file.path(HOME_TORSTAN, 'cmdstan')))

hmctdmrest <- function(drug=NULL, 
                    data=NULL,
                    prior=NULL, 
                    stan_model_option=NULL,
                    stan_sample_option=NULL
){
  object <- list(status="init")
  
  object <- object %>% 
              preprocess_valid_container( 
                                          drug, 
                                          data, 
                                          prior, 
                                          stan_model_option, 
                                          stan_sample_option)

  stan_result <- hmctdmr_stan_container(object)
  mrg.mod <- object %>% hmctdmr_mrgsolve
  

  object$stan_result <- stan_result
  object <- hmctdmr_postprocess(object)

  class(object) <- "hmctdm"

  return(object)

}
