hmctdmr_mrgsolve <- function(hmctdm, ...){
  
  mrg.mod <- get_mrgmod(hmctdm$drug)
	
}
#' get_mrgmod() Function
#'
#' This function adds a and b
#' @param drug
#' @keywords hmctdmr
#' @export
#' @examples
#' get_mrgmod()
get_mrgmod <- function(drug){
  mrg.mod <- mrgsolve::mread(drug, mrglib())
}
