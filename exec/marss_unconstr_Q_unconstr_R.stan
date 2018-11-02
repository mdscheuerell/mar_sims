data {
  // scalars
  int<lower=0> n_year;           // # years
  int<lower=0> n_spp;            // # species
  int<lower=0> n_off;            // # non-zero off-diagonals
  int<lower=1> n_q;              // # unique proc SD
  int<lower=1> n_r;              // # unique obs SD
  int<lower=1> n_obs;            // # observed ts
  // vectors
  int id_q[n_spp];                 // IDs for proc SD
  int id_r[n_obs];                 // IDs for obs SD
  // matrices
  row_vector[n_obs] yy[n_year];     // data
  int<lower=0> rc_off[n_off,2];  // indices of non-zero off-diags
  matrix<lower=0,upper=1>[n_spp,n_obs] Zmat;  // Z matrix
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
  cholesky_factor_corr[n_spp] L_corr_proc;  // chol factor
  vector<lower=0>[n_r] SD_obs;            // obs SD
  cholesky_factor_corr[n_obs] L_corr_obs;  // chol factor
  matrix[n_year,n_spp] xx;                // states
  real ymiss[n_na];
}
transformed parameters {
  // SD
  vector[n_spp] sig_proc;
  vector[n_obs] sig_obs;
  // B matrix
  matrix[n_spp,n_spp] Bmat;
  row_vector[n_obs] yymiss[n_year];     // data
  // B diagonal
  Bmat = diag_matrix(Bdiag);
  // B off-diagonals
  for(i in 1:n_off) {
    Bmat[rc_off[i,1],rc_off[i,2]] = Boffd[i];
  }
  // proc SD
  for(i in 1:n_spp) {
    sig_proc[i] = SD_proc[id_q[i]];
  }
  // obs SD
  for(i in 1:n_obs) {
    sig_obs[i] = SD_obs[id_r[i]];
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
  matrix[n_year,n_obs] xZ_tmp;
  // priors
  // intial states
  // X0 ~ normal(0,1);
  // process SD's
  SD_proc ~ cauchy(0, 3);
  L_corr_proc ~ lkj_corr_cholesky(1);
  // obs SD
  SD_obs ~ cauchy(0, 3);
  L_corr_obs ~ lkj_corr_cholesky(1);
  // B matrix
  // diagonal
  Bdiag ~ beta(1.05,1.05);
  // Bdiag ~ normal(0.5,10);
  // off-diagonals
  Boffd ~ normal(0,10);
  // likelihood
  xZ_tmp = xx * Zmat;
  for(t in 1:n_year) {
    if (t > 1) {
      xx[t] ~ multi_normal_cholesky(xx[t-1] * Bmat, diag_pre_multiply(sig_proc, L_corr_proc));
    }
    yymiss[t] ~ multi_normal_cholesky(xZ_tmp[t], diag_pre_multiply(sig_obs, L_corr_obs));
  }
}
