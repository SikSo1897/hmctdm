functions {
  vector TwoCptIV_Model(real t, vector cmt, real[] parms, real[] x_r, int[] x_i){
    
    real CL  = parms[1];
    real VC  = parms[2];
    real k12 = parms[3];
    real k21 = parms[4];
    
    real ke     = CL/VC;

    vector[2] dxdt_cmt;

    dxdt_cmt[1] = k21 * cmt[2] - (ke + k12) * cmt[1];
    dxdt_cmt[2] = k12 * cmt[1] - k21 * cmt[2];

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
  
  real<lower=0> CrCL[nt];
  real<lower=0> LBW[nt];

  real Prior_CL_NR;
  real Prior_VC_NR;
  real Prior_CL_SLOPE;
  real Prior_k12;
  real Prior_k21;

  real Prior_CL_NR_omega;
  real Prior_VC_NR_omega;
  real Prior_CL_SLOPE_omega;
  real Prior_k12_omega;
  real Prior_k21_omega;

}

transformed data{
  
  int nTheta = 4;
  int nCmt = 2;

}

parameters {
  real<lower=0> CL_NR;
  real<lower=0> VC_NR;
  real<lower=0> CL_SLOPE;
  real<lower=0> k12;
  real<lower=0> k21;
}

transformed parameters {
  real theta[nTheta];

  row_vector[nt] cHat;
  row_vector[nObs] cHatObs;

  real<lower=0> sigmaEPS[nObs];
  
  real CL = ( CL_NR *  max(LBW) + CL_SLOPE * max(CrCL) ) * 60 / 1000;
  real VC = VC_NR *  max(LBW);
  real VP = k12/k21 * VC;

  theta[1] = CL;
  theta[2] = VC;
  theta[3] = k12;
  theta[4] = k21;

  matrix[nCmt, nt] x = pmx_solve_rk45(TwoCptIV_Model, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  cHat = x[1, ] ./ VC;

  for (i in 1:nObs) {
    cHatObs[i] = cHat[iObs[i]];
    sigmaEPS[i] = sigma(cHatObs[i]);
  }
}

model {

  CL_NR ~ lognormal(log(Prior_CL_NR), Prior_CL_NR_omega);
  VC_NR ~ lognormal(log(Prior_VC_NR), Prior_VC_NR_omega);
  CL_SLOPE ~ lognormal(log(Prior_CL_SLOPE), Prior_CL_SLOPE_omega);
  k12 ~ lognormal(log(Prior_k12), Prior_k12_omega);
  k21 ~ lognormal(log(Prior_k21), Prior_k21_omega);

  for ( i in 1:nObs){
    cObs[i] ~ normal(cHatObs[i], sigmaEPS[i]);
  }
  
}