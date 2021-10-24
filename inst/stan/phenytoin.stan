functions {

  vector ODE_NonLinear(real t, vector cmt, real[] parms, real[] x_r, int[] x_i){
    
    real Vmax = parms[1];
    real Km = parms[2];
    real VD = parms[3];
    real CL = parms[4];

    vector[1] dxdt_cmt;

    real CP = cmt[1]/VD;
    
    dxdt_cmt[1] = -( CL + (Vmax/24) / (Km + CP) ) * CP;

    return dxdt_cmt;
  }
  real cv_to_sd(real cv) {
    return sqrt(log(cv^2+1));
  }

  real sigma(real conc){
    return (conc * 0.1 + 1.0);
  }

}

data {

  int<lower=1> nt;                        // number of events
  int<lower=1> nObs;
  int<lower=1> iObs[nObs];
  vector<lower=0>[nObs] cObs;
  
  // torstan parameters (nonmem data)
  int<lower=1> cmt[nt];

  int evid[nt];
  int addl[nt];
  int ss[nt];
  real amt[nt];
  real time[nt];
  real rate[nt];
  real ii[nt];

  real CLCR[nt];
  real LBW[nt];  
  real WT[nt];

}

transformed data{
  int nCmt = 1;
  real CL_SLOPE=0.01;
}
parameters{
  real ln_Vmax;
  real ln_Km;
  real ln_VD_NR;

}
transformed parameters  {

  real parms[4];

  real Vmax = exp(ln_Vmax);
  real Km = exp(ln_Km);
  real VD_NR = exp(ln_VD_NR);
  
  real CL = CL_SLOPE * CLCR[1] * 60 / 1000;
  real VD = VD_NR * LBW[1];

  row_vector[nt] cHat;
  row_vector[nObs] cHatObs;

  real<lower=0> sigmaEPS[nObs];
  
  parms[1] = Vmax * pow(WT[1]/70, 0.6);
  parms[2] = Km;
  parms[3] = VD;
  parms[4] = CL;

  matrix[nCmt, nt] x = pmx_solve_rk45(ODE_NonLinear, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, parms, 1e-5, 1e-8, 1e5);   
  
  cHat = x[1, ] ./ VD;

  for (i in 1:nObs) {
    cHatObs[i] = cHat[iObs[i]];
    sigmaEPS[i] = sigma(cHatObs[i]);
  }
}

model {

   /* ... declarations ... statements ... */

  //  VD_NR ~ lognormal(log(0.8), cv_to_sd(0.2));
  //  Km ~ lognormal(log(5.0), cv_to_sd(0.5));
  //  Vmax ~ lognormal(log(500), cv_to_sd(0.3));

   ln_VD_NR ~ normal(log(0.8), cv_to_sd(0.2));
   ln_Km ~ normal(log(5.0), cv_to_sd(0.2));
   ln_Vmax ~ normal(log(500), cv_to_sd(0.2));

   for (i in 1:nObs){
    cObs[i] ~ normal(cHatObs[i], sigmaEPS[i]);
   }

}