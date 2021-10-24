
preprocess_valid_container <- function(object, drug, data, prior, stan_model_option, stan_sample_option){
	
	object %>% 
		preprocess_valid_drug(drug) %>%
		preprocess_valid_data(data) %>%
		preprocess_valid_env(stan_model_option, stan_sample_option)
		# preprocess_valid_prior(prior) 

	# preprocess_valid_data(hmctdm, data)
	# preprocess_valid_prior(hmctdm, prior)
	# preprocess_valid_stan_model_option(hmctdm, stan_model_option)
	# preprocess_valid_stan_sample_option(hmctdm, stan_sample_option)

}

preprocess_valid_drug <- function(object, drug){
	drugs <- c("amikacin", "vancomycin", "theophiline", "phenytoin")

	if (drug %in% drugs) {
		object$drug <- drug
	} else {
		object$status <- "abort"
		stop("Unavailable Drug")
	}

	return(object)
}

preprocess_valid_data <- function(object, data){

	if (is.null(data)){

		sample.data <- load_sample_data(drug)

		message("data format")
		# print(sample.data)
		# stop("data is null")

		sample.data <- read.csv("/home/sikso/work/00_Private/hmctdm/inst/sample/sample.amikacin.csv")

		object$status <- "Sample"
		data <- sample.data

		# stop("Data is empty")
	}
	
	# print(data)
	`%notin%` <- Negate(`%in%`)

	data <- data %>%
		mutate(LBW=ifelse("LBW" %notin% colnames(data) ,
			 Peck_LBW(AGE, SEX, WT, HT), LBW), .after=WT) %>%
		mutate(CrCL=ifelse("CrCL" %notin% colnames(data) ,
			 Cockroft_Gault_CrCl(AGE=AGE, SEX=SEX, LBW=LBW, SCR=SCR), LBW), .after=SCR) 

	object$data <- data

	return (object)
}

preprocess_valid_prior <- function(object, drug, prior) {
	  stanD <- switch(drug,
                "amikacin" = function(){
			message("STAND::AMI")
			list( CL_NR=0.0417,
				CL_SLOPE=0.815,
				VD_NR = 0.27,
				VD_SLOPE=0.815) },

                "vancomycin" = function(){
			message("STAND::vancomycin")
			list( CL_NR = 0.05,
				VC_NR = 0.5, 
				CL_SLOPE=0.75,
				k12=1.12,
				k21=0.48 )},
                "theophiline" = function(){
			list( VC_NR=0.5, CL_NR=40)
		},
                "phenytoin" = file.path(stanLib(), "phenytoin.stan"),
                NULL
              )
list( CL_NR=0.0417,
                        CL_SLOPE=0.815,
                        VD_NR = 0.27,
                        VD_SLOPE=0.815 )
}

preprocess_valid_env <- function(object, stan_model_option, stan_sample_option){

	`%notin%` <- Negate(`%in%`)

	if ( "dir" %notin% names(stan_model_option) ){
		stan_model_option$dir <- file.path(getwd(), "stanc")
		if(!(dir.exists(file.path(getwd(), "stanc")))) dir.create(file.path(getwd(), "stanc"))
	}
	if ( "compile" %notin% names(stan_model_option) ){
		stan_model_option$compile <- TRUE
	}

	if ( "chains" %notin% names(stan_sample_option) ){
		stan_sample_option$chains <- 4
	}
	if ( "thin" %notin% names(stan_sample_option) ){
		stan_sample_option$thin <- 1
	}
	if ( "iter_warmup" %notin% names(stan_sample_option) ){
		stan_sample_option$iter_warmup <- 2500 * stan_sample_option$thin
	}
	if ( "iter_sampling" %notin% names(stan_sample_option) ){
		stan_sample_option$iter_sampling <- 2500 * stan_sample_option$thin 
	}
	if ( "output_dir" %notin% names(stan_sample_option) ){
		stan_sample_option$output_dir <- file.path(getwd(), "stan_output")
		if(!dir.exists(file.path(getwd(), "stan_output"))) dir.create(file.path(getwd(), "stan_output"))
	}
	if ( "adapt_delta" %notin% names(stan_sample_option) ){
		stan_sample_option$adapt_delta <- 0.95
	}

	object$stan_model_option <- stan_model_option
	object$stan_sample_option <- stan_sample_option

	return(object)
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
