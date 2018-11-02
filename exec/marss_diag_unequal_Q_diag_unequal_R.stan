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
  int<lower=0> n_pos; // number of non-missing observations
  int<lower=0> row_indx_pos[n_pos];
  int<lower=0> col_indx_pos[n_pos];
  int<lower=0> n_na; // number of missing observations
  int<lower=0> row_indx_na[n_na];
  int<lower=0> col_indx_na[n_na];
}
parameters {
  vector<lower=-1,upper=1>[n_off] Boffd;  // off-diags of B
  vector<lower=0,upper=1>[n_spp] Bdiag;   // diag of B
  vector<lower=0>[n_q] SD_proc;           // proc SD
  vector<lower=0>[n_r] SD_obs;            // obs SD
  vector[n_spp] X0;                       // initial states
  matrix[n_spp,n_year] xx;                // states
  real ymiss[n_na];
}
transformed parameters {
  // expectation of states
  matrix[n_spp,n_year] Ex;
  // B matrix
  matrix[n_spp,n_spp] Bmat;
  matrix[n_obs,n_year] yymiss;
  // diagonal
  Bmat = diag_matrix(0.5 + Bdiag);
  // off-diagonals
  for(i in 1:n_off) {
    Bmat[rc_off[i,1],rc_off[i,2]] = Boffd[i];
  }
  // expectations
  for(t in 1:n_year) {
    if (t < 2)
      Ex[,t] = Bmat * X0;
    else
      Ex[,t] = Bmat * xx[,t-1];
  }
  // Deal with missing and non-missing values separately
  for(i in 1:n_pos) {
    yymiss[row_indx_pos[i], col_indx_pos[i]] = yy[row_indx_pos[i], col_indx_pos[i]];
  }
  // Include missing observations
  if(n_na > 0) {
    for(i in 1:n_na) {
      yymiss[row_indx_na[i], col_indx_na[i]] = ymiss[i];
    }
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
  // off-diagonals
  Boffd ~ normal(0,10);
  // likelihood
  for(i in 1:n_spp) {
    row(xx,i) ~ normal(row(Ex,i), SD_proc[id_q[i]]);
  }
  for(i in 1:n_obs) {
    row(yymiss,i) ~ normal(row(Zmat * xx, i), SD_obs[id_r[i]]);
  }
}
