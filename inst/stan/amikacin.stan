functions {
  vector OneCptIV_Model(real t, vector cmt, real[] parms, real[] x_r, int[] x_i){
    
    real CL = parms[1];
    real VD = parms[2];
    real ke     = CL/VD;

    vector[1] dxdt_cmt;
    
    dxdt_cmt[1] = -ke * cmt[1];

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
  real Prior_VD_NR; 
  real Prior_CL_SLOPE;
 
  real Prior_CL_NR_omega;
  real Prior_VD_NR_omega;
  real Prior_CL_SLOPE_omega;
  
}

transformed data{
  int nTheta = 2;
  int nCmt = 1;
}

parameters {

  real CL_NR;
  real VD_NR;
  real CL_SLOPE;

}

transformed parameters {
  real theta[nTheta];

  row_vector[nt] cHat;
  row_vector[nObs] cHatObs;
  real<lower=0> sigmaEPS[nObs];
  
  real CL = (CL_NR *  max(LBW)  + CL_SLOPE * max(CrCL) ) *  60 / 1000;
  real VD = (VD_NR *  max(LBW));
  
  theta[1] = CL;
  theta[2] = VD;
    
  matrix[nCmt, nt] x = pmx_solve_rk45(OneCptIV_Model, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  cHat = x[1, ] ./ VD;

  for (i in 1:nObs) {
    cHatObs[i] = cHat[iObs[i]];
    sigmaEPS[i] = sigma(cHatObs[i]);
  }
}

model {

  CL_NR ~ lognormal(log(Prior_CL_NR), Prior_CL_NR_omega);
  VD_NR ~ lognormal(log(Prior_VD_NR), Prior_VD_NR_omega);
  CL_SLOPE ~ lognormal(log(Prior_CL_SLOPE), Prior_CL_SLOPE_omega);

  for ( i in 1:nObs){
    cObs[i] ~ normal(cHatObs[i], sigmaEPS[i]);
  }
}