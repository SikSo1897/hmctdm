functions {
  vector OneCptOral_Model(real t, vector cmt, real[] parms, real[] x_r, int[] x_i){
    
    real CL = parms[1];
    real VC = parms[2];
    real ka = parms[3];
    
    real ke     = CL/VC;

    vector[2] dxdt_cmt;
    
    dxdt_cmt[1] = -ka * cmt[1];
    dxdt_cmt[2] = ka * cmt[1] - ke * cmt[2];

    return dxdt_cmt;
  }

  real cv_to_sd(real cv) {
    return sqrt(log(cv^2+1));
  }

  real sigma(real conc){
    return (conc * 0.15 + 0.25);
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
  
  real<lower=0> CLCR[nt];
  real<lower=0> LBW[nt];
  
}

transformed data{
  
  // @ Prior
  real TV_VC_NR = 0.5; 
  real TV_CL_NR = 40;

  real OMEGA_VC_NR = cv_to_sd(0.2);
  real OMEGA_CL_NR = cv_to_sd(0.15);
  
  real KA = 0.27;

  int nTheta = 3;
  int nCmt = 2;

}

parameters {
  real<lower=0> CL_NR;
  real<lower=0> VC_NR;
}

transformed parameters {
  real theta[nTheta];

  row_vector[nt] cHat;
  row_vector[nObs] cHatObs;
  real<lower=0> sigmaEPS[nObs];
  
  real CL = CL_NR *  max(LBW) / 1000;
  real VC = VC_NR *  max(LBW);
  
  theta[1] = CL;
  theta[2] = VC;
  theta[3] = KA;

  matrix[nCmt, nt] x = pmx_solve_rk45(OneCptOral_Model, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  cHat = x[2, ] ./ VC;

  for (i in 1:nObs) {
    cHatObs[i] = cHat[iObs[i]];
    sigmaEPS[i] = sigma(cHatObs[i]);
  }
}

model {

  CL_NR ~ lognormal(log(TV_CL_NR), OMEGA_CL_NR);
  VC_NR ~ lognormal(log(TV_VC_NR), OMEGA_VC_NR);

  for ( i in 1:nObs){
    cObs[i] ~ normal(cHatObs[i], sigmaEPS[i]);
  }
}