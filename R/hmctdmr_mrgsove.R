hmctdmr_mrgsolve <- function(hmctdm, ...){
  
  mrg.mod <- get_mrgmod(hmctdm$drug)
	
}

get_mrgmod <- function(drug){
  mrg.mod <- mrgsolve::mread(drug, mrglib())
}
