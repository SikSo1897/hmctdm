
preprocess_valid_container <- function(object, drug, data, prior, stan_model_option, stan_sample_option){
	
	object %>% 
		preprocess_valid_drug(drug) %>%
		preprocess_valid_data(data) %>%
		preprocess_valid_env(stan_model_option, stan_sample_option) %>%
		preprocess_valid_prior(prior) 
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

		sample.data <-get_sample_data(object$drug)

		object$status <- "Sample"
		data <- sample.data

	}
	
	`%notin%` <- Negate(`%in%`)

	data <- data %>%
		mutate(LBW=ifelse("LBW" %notin% colnames(data) ,
			 Peck_LBW(AGE, SEX, WT, HT), LBW), .after=WT) %>%
		mutate(CrCL=ifelse("CrCL" %notin% colnames(data) ,
			 Cockroft_Gault_CrCl(AGE=AGE, SEX=SEX, LBW=LBW, SCR=SCR), LBW), .after=SCR) 

	object$data <- data

	return (object)
}

preprocess_valid_prior <- function(object, prior) {
	`%notin%` <- Negate(`%in%`)

	new_prior <- NULL

	x <- get_prior(object$drug)
	
	check <- (x$name %notin% names(prior))

	for ( i in 1:length(x$name)){
		
		if (check[i]){
			new_prior <- append(new_prior, x$value[i])
		} else {
			new_prior <- append(new_prior, with(prior, get(x$name[i])))
			
		}
	}

	new_prior <- setNames(as.list(new_prior), x$stan_prior_name) 

	object$prior <- new_prior

	return(object)

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

get_sample_data <- function(drug){
	read.csv(base::file.path(sampleLib(), paste0(drug, ".csv")))
}

get_prior <- function(drug){
		x <- 	switch(drug,

		"amikacin" = list(
										name = c("CL_NR", "VD_NR", "CL_SLOPE", 
															"CL_NR_omega", "VD_NR_omega", "CL_SLOPE_omega"),

										stan_prior_name = c("Prior_CL_NR", "Prior_VD_NR", "Prior_CL_SLOPE", 
																					"Prior_CL_NR_omega", "Prior_VD_NR_omega", "Prior_CL_SLOPE_omega"),

										value = c(0.0417, 0.27, 0.815, 
															cv_to_sd(0.25), cv_to_sd(0.3), cv_to_sd(0.4))
									),
									
		"vancomycin"= list(
										name = c("CL_NR", "VC_NR", "CL_SLOPE", "k12", "k21",
														 	"CL_NR_omega", "VC_NR_omega", "CL_SLOPE_omega", "k12_omega", "k21_omega"),

										stan_prior_name = c("Prior_CL_NR", "Prior_VC_NR", "Prior_CL_SLOPE", "Prior_k12", "Prior_k21",
																					"Prior_CL_NR_omega", "Prior_VC_NR_omega", "Prior_CL_SLOPE_omega", "Prior_k12_omega", "Prior_k21_omega"),

										value = c(0.05, 0.21, 0.75, 1.12, 0.48,
															cv_to_sd(0.15), cv_to_sd(0.2), cv_to_sd(0.33), cv_to_sd(0.25), cv_to_sd(0.25))
									),

		"theophiline"= list(
										name = c("ka", "CL_NR", "VC_NR",
																"CL_NR_omega", "VC_NR_omega"),

										stan_prior_name = c("Prior_ka", "Prior_CL_NR", "Prior_VC_NR",
																					"Prior_CL_NR_omega", "Prior_VC_NR_omega"),

										value = c(0.27, 40, 0.5,
																cv_to_sd(0.15), cv_to_sd(0.2))
									),

		"phenytoin"= list(
										name = c("Prior_Vmax", "Prior_Km", "Prior_VD_NR", 
																"Prior_Vmax_omega", "Prior_Km_omega", "Prior_VD_NR_omega"),

										stan_prior_name = c("Prior_Vmax", "Prior_Km", "Prior_VD_NR", 
																				"Prior_Vmax_omega", "Prior_Km_omega", "Prior_VD_NR_omega"),

										value = c(500, 5.0, 0.8,
															cv_to_sd(0.2), cv_to_sd(0.2), cv_to_sd(0.2))
									))
}

get_default_prior <- function(drug){
	
	default_prior <- NULL

	x <- get_prior(drug)
	

	for ( i in 1:length(x$name)){
			default_prior <- append(default_prior, x$value[i])
	}

	default_prior <- setNames(as.list(default_prior), x$stan_prior_name) 
	return(default_prior)

}