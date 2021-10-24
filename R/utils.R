#' @noRd
mrglib <- function(){
  system.file("mrgsolve", package = "hmctdmr")

}

#' @noRd
stanLib <- function(){
  system.file("stan", package="hmctdmr")

}

#' @noRd
sampleLib <- function(){
  system.file("sample", package="hmctdmr")

}


Peck_LBW <- function(age, sex, wt, ht) { 
	ht  <- ht/2.54 
	LBW <- ifelse(sex == 0, 
		ifelse(age <= 18,
			-59.6035 + (5.2878 * ht) - (0.123939 * ht ** 2) + (0.00128936 * ht ** 3),
			-130.736 + (4.064 *ht)),
		ifelse(age <= 18,
			-77.55796 + (6.93728 * ht) - (0.171703 * ht ** 2) + (0.001726 * ht ** 3),
			-111.621 + (3.636 * ht))
	)/(2.204623)
	ifelse(LBW>=wt, wt, LBW) %>% return
}
Cockroft_Gault_CrCl <- function(AGE, SEX, LBW, SCR, crcl_lower=10.3, crcl_upper=127.1){
	# return(1)
	crcl <- (140-AGE) * LBW / 72 * SCR
	crcl <- ifelse(SEX==0, crcl, crcl * 0.85) 
	crcl <- ifelse(crcl >= crcl_upper, crcl_upper, crcl)
	crcl <- ifelse(crcl <= crcl_lower, crcl_lower, crcl)
	return(crcl)
}

load_sample_data <- function(drug){
  dplyr::tibble(ID=1, amt=1000,)
}


cv_to_sd <- function(cv) sqrt(log(cv^2 + 1))