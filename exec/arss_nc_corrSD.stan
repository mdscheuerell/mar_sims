data {
  // length of ts
  int<lower=0> n_year;
  // data
  vector[n_year] yy;
  // mean of log-normal
  vector[2] mean_SD;
}
parameters {
  // proc SD
  // real<lower=0> SD_proc;
  // obs SD
  // real<lower=0> SD_obs;
  // MVN for log(SD)
  cholesky_factor_corr[2] Lcorr;
  real<lower=0> sigma_SD;
  vector[2] SD_vec;
  // AR(1) coef
  real<lower=0, upper=1> phi;
  // initial state
  real init_state;
  // unit normals
  vector[n_year] zz;
}
transformed parameters {
  // proc SD
  real<lower=0> SD_proc;
  // obs SD
  real<lower=0> SD_obs;
  vector<lower=0>[2] sigma_vec;
  // states
  vector[n_year] xx;
  SD_proc = SD_vec[1];
  SD_obs = SD_vec[2];
  sigma_vec[1] = sigma_SD;
  sigma_vec[2] = sigma_SD;
  // initial state
  xx[1] = phi * init_state + SD_proc * zz[1];
  // remaining states
  for(t in 2:n_year) {
    xx[t] = phi * xx[t-1] + SD_proc * zz[t];
  }
}
model {
  // PRIORS
  // process SD
  // SD_proc ~ normal(0, 1);
  // obs SD
  // SD_obs ~ normal(0, 1);
  SD_vec ~ multi_normal_cholesky(mean_SD, diag_pre_multiply(sigma_vec, Lcorr));
  // mean_SD ~ normal(0, 1);
  sigma_SD ~ normal(0, 1);
  Lcorr ~ lkj_corr_cholesky(1);
  // unit normals
  zz ~ std_normal();
  // AR(1) coef
  phi ~ normal(0.5, 2);
  // initial state
  init_state ~ normal(0, SD_proc/(1 - phi^2));
  // LIKELIHOOD
  yy ~ normal(xx, SD_obs);
}
generated quantities {
  matrix[2,2] Sigma;
  Sigma = diag_pre_multiply(sigma_vec, Lcorr);
}
