### Installation

You can install development version from github by executing the following code in R console.

``` r
install.packages("devtools")
devtools::install_github("sikso1897/hmctdm")
```

hmctdm relies on 
[torsten](#)
...

cmdstanr path

``` r
HOME_TORSTAN <- "path/to/torsten"
suppressMessages(set_cmdstan_path(file.path(HOME_TORSTAN, 'cmdstan')))
```

### Example
``` r
library(cmdstanr)
library(hmctdm)
```

#### 1) set cmdstan path
``` r
HOME_TORSTAN <- "path/to/torsten"
suppressMessages(set_cmdstan_path(file.path(HOME_TORSTAN, 'cmdstan')))
```
#### 2) bring data set
```r

```

#### 3) estimate
```r
hmctdm <- hmctdmr::hmctdmrest(drug="amikacin")
print(hmctdm)
```

``` r
# search stan_file ..
# find stan_file : /home/sikso/work/00_Private/hmctdm/inst/stan/amikacin.stan
# Compiling Stan program...
# Running MCMC with 4 sequential chains...
# 
# Chain 1 finished in 0.6 seconds.
# Chain 2 finished in 0.5 seconds.
# Chain 3 finished in 0.6 seconds.
# Chain 4 finished in 0.6 seconds.
# 
# All 4 chains finished successfully.
# Mean chain execution time: 0.6 seconds.
# Total execution time: 3.1 seconds.
# Building amikacin ... done.
#   ID time evid amt cmt ss ii addl rate SEX AGE       WT      LBW       HT      SCR
# 1  1    0    1 500   1  0  0    0 1000   1  79 43.97928 43.97928 155.7253 1.072065
# 2  1    1    0   0   1  0  0    0    0   1  79 43.97928 43.97928 155.7253 1.072065
#       CrCL      DV
# 1 33.95356  0.0000
# 2 33.95356 28.5206
# 
# Status:  Sample 
# Drug:  Amikacin 
# 
# 
# Stan sample option: 
# # A tibble: 1 × 3
#   option            dir                                           compile
#   <chr>             <chr>                                         <chr>  
# 1 stan_model_option /home/sikso/work/00_Private/hmctdm/test/stanc TRUE   
# 
# Stan model option: 
# # A tibble: 1 × 7
#   option             chains thin  iter_warmup iter_sampling output_dir   adapt_delta
#   <chr>              <chr>  <chr> <chr>       <chr>         <chr>        <chr>      
# 1 stan_sample_option 4      1     2500        2500          /home/sikso… 0.95       
# 
# Estimated PK Parameters: 
# # A tibble: 5 × 10
#   variable    mean  median     sd    mad      q5     q95  rhat ess_bulk ess_tail
#   <chr>      <dbl>   <dbl>  <dbl>  <dbl>   <dbl>   <dbl> <dbl>    <dbl>    <dbl>
# 1 CL_NR     0.0430  0.0416 0.0110 0.0102  0.0276  0.0627  1.00    4780.    3047.
# 2 VD_NR     0.337   0.335  0.0558 0.0571  0.246   0.430   1.00    5766.    4907.
# 3 CL_SLOPE  0.932   0.856  0.389  0.340   0.443   1.68    1.00    4380.    3303.
# 4 CL        2.01    1.86   0.794  0.693   1.02    3.54    1.00    4417.    3270.
# 5 VD       14.8    14.7    2.46   2.51   10.8    18.9     1.00    5766.    4907.
# 
# Estimated Concentration: 
#   ID time evid amt cmt ss ii addl rate SEX AGE       WT      LBW       HT      SCR
# 1  1    0    1 500   1  0  0    0 1000   1  79 43.97928 43.97928 155.7253 1.072065
# 2  1    1    0   0   1  0  0    0    0   1  79 43.97928 43.97928 155.7253 1.072065
#       CrCL      DV     cHat
# 1 33.95356  0.0000  0.00000
# 2 33.95356 28.5206 30.52245
# 
# Individual Data: 
#   ID SEX AGE       WT      LBW       HT      SCR     CrCL     CL_NR     VD_NR
# 1  1   1  79 43.97928 43.97928 155.7253 1.072065 33.95356 0.0415959 0.3349125
#    CL_SLOPE      CL      VD
# 1 0.8558655 1.85598 14.7292

```
