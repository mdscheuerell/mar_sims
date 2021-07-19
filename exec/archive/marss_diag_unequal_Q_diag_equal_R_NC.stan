data {
  // scalars
  int<lower=0> n_year;           // # years
  int<lower=0> n_spp;            // # species
  int<lower=0> n_off;            // # non-zero off-diagonals
  int<lower=1> n_q;              // # unique proc SD
  int<lower=1> n_obs;            // # observed ts
  // vectors
  int id_q[n_spp];                 // IDs for proc SD
  // matrices
  int<lower=0> rc_off[n_off,2];  // indices of non-zero off-diags
  matrix<lower=0,upper=1>[n_obs,n_spp] Zmat;
  int<lower=0> n_pos; // number of non-missing observations
  int<lower=0> row_indx_pos[n_pos];
  int<lower=0> col_indx_pos[n_pos];
  int<lower=0> n_na; // number of missing observations
  int<lower=0> row_indx_na[n_na];
  int<lower=0> col_indx_na[n_na];
  real yy[n_pos];       // data
}
parameters {
  real<lower=0> SD_obs;
  vector<lower=-1,upper=1>[n_off] Boffd;  // off-diags of B
  vector<lower=0,upper=1>[n_spp] Bdiag;   // diag of B
  // vector<lower=0>[n_q] SD_proc;           // proc SD
  real<lower=0> SD_proc;           // proc SD
  vector[4] init_state;           // proc SD
  matrix[n_spp,n_year] Z_proc;
  // real<lower=0> var_SD;
  // real<upper=0> cov_SD;
  // cholesky_factor_corr[2] Lcorr;
  // real<lower=0> sigma;
  // vector[2] SD_vec;
  // matrix[n_spp,n_year] xx;       // states
  vector[n_na] ymiss;  // missing values
}
transformed parameters {
  // obs SD
  // real<lower=0> SD_obs;
  // real<lower=0> SD_proc;
  // vector<lower=0>[2] sigma_vec;
  // vector<lower=0>[2] mean_SD;
  // cov_matrix[2] sigma_SD;
  // cov matrix
  // cov_matrix[n_spp] QQ;
  // B matrix
  matrix[n_spp,n_spp] Bmat;
  matrix[n_obs,n_year] yymiss;
  matrix[n_spp,n_year] xx;       // states
  // diagonal
  Bmat = diag_matrix(Bdiag);
  // off-diagonals
  for(i in 1:n_off) {
    Bmat[rc_off[i,1],rc_off[i,2]] = Boffd[i];
  }
  // cov matrix
  // for(i in 1:n_spp) {
  // 	QQ[i,i] = SD_proc[id_q[i]]^2;
  // }
  // for(i in 1:(n_spp-1)) {
  // 	for(j in (i+1):n_spp) {
  // 	  QQ[i,j] = 0;
  // 	  QQ[j,i] = 0;
  // 	}
  // }
  // Deal with missing and non-missing values separately
  for(i in 1:n_pos) {
    yymiss[row_indx_pos[i], col_indx_pos[i]] = yy[i];
  }
  // Include missing observations
  if(n_na > 0) {
    for(i in 1:n_na) {
      yymiss[row_indx_na[i], col_indx_na[i]] = ymiss[i];
    }
  }
  // SD_obs = exp(SD_vec[1]);
  // SD_proc = exp(SD_vec[2]);
  // mean_SD[1] = 0;
  // mean_SD[2] = 0;
  // sigma_vec[1] = sigma;
  // sigma_vec[2] = sigma;
  // sigma_SD[1,1] = var_SD;
  // sigma_SD[1,2] = cov_SD;
  // sigma_SD[2,1] = cov_SD;
  // sigma_SD[2,2] = var_SD;
  xx[,1] = init_state;
  for(t in 2:n_year) {
    xx[,t] = Bmat * xx[,t-1] + SD_proc * Z_proc[,t];
  }
}
model {
  // PRIORS
  // initial state
  // xx[,1] ~ normal(0,5);
  init_state ~ normal(0,5);
  // process SD's
  // SD_proc ~ normal(0,1);
  to_vector(Z_proc) ~ normal(0,1);
  // obs SD
  SD_obs ~ normal(0,1);
  // SD_vec ~ multi_normal_cholesky(mean_SD, diag_pre_multiply(sigma_vec, Lcorr));
  // sigma ~ normal(0,1);
  // Lcorr ~ lkj_corr_cholesky(0.8);
  // B matrix
  // diagonal
  // Bdiag ~ beta(1.2,1.2);
  Bdiag ~ normal(0.5,1);
  // off-diagonals
  Boffd ~ normal(0,1);
  // missing obs
  ymiss ~ normal(0,5);
  // LIKELIHOOD
  // for(t in 2:n_year) {
  //   xx[,t] ~ normal(Bmat * xx[,t-1], SD_proc);
  //   // xx[,t] = Bmat * xx[,t-1] + SD_proc * Z_proc;
  //   // col(xx,t) ~ multi_normal(Bmat * col(xx,t-1), QQ);
  // }
  // to_vector(yymiss) ~ normal(to_vector(xx), SD_obs);
  to_vector(yymiss) ~ normal(to_vector(xx), SD_obs);
}

