data {
  // scalars
  int<lower=0> n_year;           // # years
  int<lower=0> n_spp;            // # species
  int<lower=1> n_q;              // # unique proc SD
  int<lower=1> n_r;              // # unique obs SD
  int<lower=1> n_obs;            // # observed ts
  // vectors
  int id_q[n_spp];               // IDs for proc SD
  int id_r[n_obs];               // IDs for obs SD
  // matrices
  matrix[n_obs,n_year] yy;       // data
  matrix<lower=0,upper=1>[n_obs,n_spp] Zmat;
}
parameters {
  vector<lower=0,upper=1>[n_spp] Bdiag;   // diag of B
  vector<lower=0>[n_q] SD_proc;           // proc SD
  // vector<lower=0>[n_r] SD_obs;            // obs SD
  real<lower=0> SD_obs;            // obs SD
  vector[n_spp] X0;                       // initial states
  matrix[n_spp,n_year] xx;                // states  
}
transformed parameters {
  // expectation of states
  matrix[n_spp,n_year] Ex;         
  // expectations
  for(t in 1:n_year) {
    if (t < 2)
      Ex[,t] = diag_matrix(Bdiag) * X0;
    else
      Ex[,t] = diag_matrix(Bdiag) * xx[,t-1];
  }
}
model {
  // priors
  // intial states
  X0 ~ normal(0,5);
  // process SD's
  SD_proc ~ normal(0, 1);
  // obs SD
  SD_obs ~ normal(0, 1);
  // B matrix
  // diagonal
  // Bdiag ~ beta(1.05,1.05);
  Bdiag ~ normal(0.5,10);
  // likelihood
  for(i in 1:n_spp) {
    row(xx,i) ~ normal(row(Ex,i), SD_proc[id_q[i]]);
  }
  for(i in 1:n_obs) {
    row(yy,i) ~ normal(row(Zmat * xx, i), SD_obs);
  }
}
