data {
  // scalars
  int<lower=0> n_year;           // # years
  // vectors
  real yy[n_year];       // data
}
parameters {
  real<lower=0> SD_obs;
  real<lower=0> SD_proc;           // proc SD
  real<lower=0.01, upper=0.99> phi;           // AR(1) coef
  real xx[n_year];       // states
}
// transformed parameters {
// }
model {
  // PRIORS
  // process SD's
  // SD_proc ~ normal(0, 1);
  SD_proc ~ uniform(0.1, 0.5);
  // obs SD
  // SD_obs ~ normal(0, 1);
  SD_obs ~ uniform(0.1, 0.5);
  // initial state
  xx[1] ~ normal(0, SD_proc/(1 - phi^2));
  // AR(1) coef
  phi ~ normal(0.5, 5);
  // LIKELIHOOD
  for(t in 2:n_year) {
    xx[t] ~ normal(phi * xx[t-1], SD_proc);
  }
  yy ~ normal(xx, SD_obs);
}

