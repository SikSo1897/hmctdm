
# print <- function(hmctdm, ...) UseMethod("pprint")
print.hmctdm <- function(hmctdm, ...){

  cat("\n")
  cat("Status: ", hmctdm$status, "\n")
  cat("Drug: ", stringr::str_to_title(hmctdm$drug), "\n\n")

  cat("\n")
  cat("Stan sample option: \n")
  print(hmctdm$stan_model_option)

  cat("\n")
  cat("Stan model option: \n")
  print(hmctdm$stan_sample_option)
  
  cat("\n")
  cat("Estimated PK Parameters: \n")
  print(hmctdm$param)

  cat("\n")
  cat("Estimated Concentration: \n")
  print(hmctdm$est)

  cat("\n")
  cat("Individual Data: \n")
  print(hmctdm$idata)

}

# plot.hmctdm <- function(object, ...){}
# hist.hmctdm <- function(object, ...){}