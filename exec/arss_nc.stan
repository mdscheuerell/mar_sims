data {
  // length of ts
  int<lower=0> n_year;
  // data
  vector[n_year] yy;
}
parameters {
  // proc SD
  real<lower=0> SD_proc;
  // obs SD
  real<lower=0> SD_obs;
  // AR(1) coef
  real<lower=0, upper=1> phi;
  // initial state
  real init_state;
  // unit normals
  vector[n_year] zz;
}
transformed parameters {
  // states
  vector[n_year] xx;
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
  SD_proc ~ normal(0, 1);
  // obs SD
  SD_obs ~ normal(0, 1);
  // unit normals
  zz ~ std_normal();
  // AR(1) coef
  phi ~ normal(0.5, 2);
  // initial state
  init_state ~ normal(0, SD_proc/(1 - phi^2));
  // LIKELIHOOD
  yy ~ normal(xx, SD_obs);
}
