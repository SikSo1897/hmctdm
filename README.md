# hmctdm 

<!-- badges: start -->
[![](https://img.shields.io/badge/Citation-NULL-blue.svg)](#)

hmctdm is open source package for bayesian estimation of PK parameters using stan in R. The main function, `hmctdmrest()` , can estimate individual PK parameters. below is a drug that can be estimate from the function:

- Amikacin (IV)
- Vancomycin (IV)
- Theophiline (Oral)
- Phenytoin (Oral)

The pharmacokinetic model was implemented from the PKS program. 

## Installation

You can install development version from github by executing the following code in R console.

``` r
install.packages("devtools")
devtools::install_github("sikso1897/hmctdm")
```

hmctdm relies on the following software & packages. 
- [Stan](https://mc-stan.org/) is an open source probabilistic programing language designed primarily to do Bayesian data analysis
- [Torsten](https://metrumresearchgroup.github.io/Torsten/) is a collection of Stan functions to facilitate analysis of pharmacometric data. 
- [Mrgsolve](https://github.com/metrumresearchgroup/mrgsolve) is model implementation and ordinary differential equation solver


Please refer to the installation guide of these programs
- [install guide of stan](https://mc-stan.org/cmdstanr/)
- [install guide of torsten](https://metrumresearchgroup.github.io/Torsten/installation/)
- [install guide of mrgsolve](https://github.com/metrumresearchgroup/mrgsolve/wiki/mrgsolve-Installation)

## Example
``` r
library(cmdstanr)
library(hmctdm)
```
#### 1) set cmdstan path
``` r
HOME_TORSTEN <- "path/to/torsten"
suppressMessages(set_cmdstan_path(file.path(HOME_TORSTEN, 'cmdstan')))
```
#### 2) bring data set
checking the sample data and bring data of type `data.frame` or `tibble` by format of sample data.
```r
sample_data <- hmctdmr::get_sample_data(drug="amikacin")
sample_data

#   ID time evid amt cmt ss ii addl rate SEX AGE       WT       HT      SCR      DV
# 1  1    0    1 500   1  0  0    0 1000   1  79 43.97928 155.7253 1.072065  0.0000
# 2  1    1    0   0   1  0  0    0    0   1  79 43.97928 155.7253 1.072065 28.5206
```

#### 3) estimate
```r
hmctdm <- hmctdmr::hmctdmrest(drug="amikacin", data="your_data_set")
print(hmctdm)
```

``` r
# Model executable is up to date!
# Running MCMC with 4 sequential chains...

# Chain 1 finished in 0.6 seconds.
# Chain 2 finished in 0.6 seconds.
# Chain 3 finished in 0.5 seconds.
# Chain 4 finished in 0.5 seconds.

# All 4 chains finished successfully.
# Mean chain execution time: 0.6 seconds.
# Total execution time: 3.0 seconds.
# Building amikacin ... done.
# Building amikacin ... (waiting) ...
# done.

# Status:  Sample 
# Drug:  Amikacin 


# Stan sample option: 
# # A tibble: 1 × 3
#   option            dir                                           compile
#   <chr>             <chr>                                         <chr>  
# 1 stan_model_option /home/sikso/work/00_Private/hmctdm/test/stanc TRUE   

# Stan model option: 
# # A tibble: 1 × 7
#   option             chains thin  iter_warmup iter_sampling output_dir                                          adapt_delta
#   <chr>              <chr>  <chr> <chr>       <chr>         <chr>                                               <chr>      
# 1 stan_sample_option 4      1     2500        2500          /home/sikso/work/00_Private/hmctdm/test/stan_output 0.95       

# Prior Information: 
# # A tibble: 6 × 2
#   parameter             value
#   <chr>                 <dbl>
# 1 Prior_CL_NR          0.0417
# 2 Prior_VD_NR          0.27  
# 3 Prior_CL_SLOPE       0.815 
# 4 Prior_CL_NR_omega    0.246 
# 5 Prior_VD_NR_omega    0.294 
# 6 Prior_CL_SLOPE_omega 0.385 

# Estimated PK Parameters: 
# # A tibble: 5 × 10
#   variable    mean  median     sd    mad      q5     q95  rhat ess_bulk ess_tail
#   <chr>      <dbl>   <dbl>  <dbl>  <dbl>   <dbl>   <dbl> <dbl>    <dbl>    <dbl>
# 1 CL_NR     0.0431  0.0419 0.0108 0.0104  0.0277  0.0626  1.00    6045.    4133.
# 2 VD_NR     0.336   0.335  0.0556 0.0554  0.247   0.429   1.00    6087.    5655.
# 3 CL_SLOPE  0.925   0.853  0.382  0.330   0.450   1.64    1.00    4353.    3045.
# 4 CL        2.00    1.85   0.778  0.673   1.03    3.45    1.00    4338.    3091.
# 5 VD       14.8    14.7    2.44   2.44   10.9    18.9     1.00    6087.    5655.

# Estimated Concentration: 
#   ID time evid amt cmt ss ii addl rate SEX AGE       WT      LBW       HT      SCR     CrCL      DV     cHat
# 1  1    0    1 500   1  0  0    0 1000   1  79 43.97928 43.97928 155.7253 1.072065 33.95356  0.0000  0.00000
# 2  1    1    0   0   1  0  0    0    0   1  79 43.97928 43.97928 155.7253 1.072065 33.95356 28.5206 30.65465

# Individual Data: 
#   ID SEX AGE       WT      LBW       HT      SCR     CrCL     CL_NR    VD_NR CL_SLOPE      CL       VD
# 1  1   1  79 43.97928 43.97928 155.7253 1.072065 33.95356 0.0418793 0.334881 0.853423 1.85247 14.72785

```
## Specific
```r
hmctdmrest <- function(drug=NULL, 
                    data=NULL,
                    prior=NULL, 
                    stan_model_option=NULL,
                    stan_sample_option=NULL){
                    ...
}
```
hmctdmrest can be given the following parameters:
#### `drug`
  - amikacin
  - vancomycin
  - theophiline
  - phenytoine
  
#### `data`

#### `prior`
It can be changed by passing the prior parameter in the list.
```r
# example

hmctdmr::get_default_prior(drug="amikacin")
# $Prior_CL_NR
# [1] 0.0417

# $Prior_VD_NR
# [1] 0.27

# $Prior_CL_SLOPE
# [1] 0.815

# $Prior_CL_NR_omega
# [1] 0.2462207

# $Prior_VD_NR_omega
# [1] 0.2935604

# $Prior_CL_SLOPE_omega
# [1] 0.3852532


hmctdm <- hmctdmr::hmctdmrest(
                    drug="amikacin", 
                    data="your_data_set",
                    prior=list(
                      CL_NR=0.5,
                      CL_NR_omega=0.3
                    ))
}
```
#### `stan_model_option` 

stan_model_option is an option for creating a new CmdStanModel object from a file containing a Stan program. More information on the options can be found [here](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html) 

#### `stan_sample_option`
stan_sample_option is option of `$sample()` method CmdStanModel objcet runs the default MCMC algorithm in CmdStan (algorithm = hmc, engine=nuts). More information on the options can be found [here](https://mc-stan.org/cmdstanr/reference/model-method-sample.html)

## Development
hmctdm is under development. Your feedback for additional feature requests or bug reporting is welcome. Contact us through the issue tracker.

