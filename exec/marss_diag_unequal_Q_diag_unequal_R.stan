data {
  // scalars
  int<lower=0> n_year;           // # years
  int<lower=0> n_spp;            // # species
  int<lower=0> n_off;            // # non-zero off-diagonals
  int<lower=1> n_q;              // # unique proc SD
  int<lower=1> n_r;              // # unique obs SD
  int<lower=1> n_obs;            // # observed ts
  // vectors
  int id_q[n_spp];               // IDs for proc SD
  int id_r[n_obs];               // IDs for obs SD
  // matrices
  matrix[n_obs,n_year] yy;       // data
  matrix<lower=0,upper=1>[n_obs,n_spp] Zmat;
  int<lower=0> rc_off[n_off,2];  // indices of non-zero off-diags
}
parameters {
  vector<lower=-1,upper=1>[n_off] Boffd;  // off-diags of B
  vector<lower=0,upper=1>[n_spp] Bdiag;   // diag of B
  vector<lower=0>[n_q] SD_proc;           // proc SD
  vector<lower=0>[n_r] SD_obs;            // obs SD
  // vector[n_spp] X0;                    // initial states
  matrix[n_spp,n_year] xx;                // states  
}
transformed parameters {
  // cov matrices
  vector[n_spp] Q;
  vector[n_obs] R;
  // B matrix
  matrix[n_spp,n_spp] Bmat;
  // diagonal
  Bmat = diag_matrix(Bdiag);
  // off-diagonals
  for(i in 1:n_off) {
    Bmat[rc_off[i,1],rc_off[i,2]] = Boffd[i];
  }
  // cov matrix Q
  for(i in 1:n_spp) {
  	Q[i] = SD_proc[id_q[i]];
  }
  // cov matrix R
  for(i in 1:n_obs) {
  	R[i] = SD_obs[id_r[i]];
  }
}
model {
  // priors
  // intial states
  // X0 ~ normal(0,1);
  // process SD's
  for(i in 1:n_q) {
  	SD_proc[i] ~ cauchy(0, 5);
  }
  // obs SD
  for(i in 1:n_r) {
  	SD_obs[i] ~ cauchy(0, 5);
  }
  // B matrix
  // diagonal
  Bdiag ~ beta(1.05,1.05);
  // off-diagonals
  Boffd ~ normal(0,10);
  // likelihood
  for(t in 2:n_year) {
    col(xx,t) ~ multi_normal(Bmat * col(xx,t-1), QQ);
    // col(yy,t) ~ multi_normal(Zmat * col(xx,t), RR);
  }
  for(n in 1:n_obs) {
    row(yy,n) ~ normal(row(Zmat * xx, n), RR);
  }
  
}
