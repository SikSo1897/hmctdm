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

#' get_recommended_dose() Function
#'
#' This function adds a and b
#' @param 
#' @keywords test
#' @export
#' @examples
#' get_recommended_dose()

get_recommended_dose <- function(
	mode, model=NULL, target=NULL,
	current_dose=NULL, current_status=NULL, ...){
	
	msg <- list()

	msg$mode <- mode
	msg$target <- target

	available_mode <- c("Cp", "b_Cp", "AUC")
	if(is.null(mode)){
		stop("Error in mode missing")
	}

	if (!(mode %in% available_mode)){
		msg <- sprintf('unavailable mode("%s") \n-- available mode list: %s', mode, toString(available_mode))
		stop(msg)
	}

	if (mode == "AUC"){

		msg$current_dose <- current_dose
		msg$current_status <- current_status

		recommended_dose <- (current_dose/current_status) * target 
	}

	if (mode == "Cp"){

		msg$current_dose <- current_dose
		msg$current_status <- current_status

		recommended_dose <- (current_dose/current_status) * target 
	}

	if (mode == "b_Cp"){
		args <- list(...)

		for (key in names(args)){
			msg[[key]] <- args[[key]]
		}

		recommended_dose <- get_recommend_dose_use_b_cp(target, ...)
	
	}
	for (key in names(msg)){
		cat("-", key, ":", msg[[key]], "\n")
	}
	cat("\n")

	cat("-","recommended_dose :", recommended_dose)
}

get_recommend_dose_use_b_cp <- function(target, Vmax, CL_R, Km, F, tau){
	recommended_dose <- (Vmax/24 + CL_R * (Km - target)) / ((Km + target) * F) * target * tau
	
	return(recommended_dose)
}

#' get_b_Cp() Function
#'
#' This function adds a and b
#' @param 
#' @keywords test
#' @export
#' @examples
#' get_b_Cp()
get_b_Cp <- function(
		dose, F, tau, CLp, 
		Vd, Vmax, Km, salt_factor=1)
{
			
	tryCatch({
		args <- list(
			dose=dose, 
			F=F, 
			tau=tau, 
			CLp=CLp, 
			Vd=Vd, 
			Vmax=Vmax, 
			Km=Km
		)
	},
	error = function(e){
		# message(e,"\n")
		stop("required arguments: {dose, F, tau, CLp, Vd, Vmax, Km} \n\n  > ", e)
	},
	warning=function(w){
		message(w,"\n")
	})
	
	r <- dose * salt_factor * F / tau
	a <- CLp/Vd
	b <- ((r-Vmax)-(CLp-Km))/Vd
	c <- Km * r / Vd

	b_Cp <- (-b + sqrt(b^2 - 4*a*c) )/(2*a)
}