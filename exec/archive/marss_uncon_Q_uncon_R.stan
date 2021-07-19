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
  vector[n_obs] yy[n_year];       // data
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
  vector[n_spp] X0;                    // initial states
  vector[n_spp] xx[n_year];                // states
  real ymiss[n_na];
}
transformed parameters {
  // temp matrices
  vector[n_obs] Ey[n_year];       // expectation of data
  vector[n_spp] Ex[n_year];       // expectation of states
  // cov matrices
  cov_matrix[n_spp] QQ;
  cov_matrix[n_obs] RR;
  // B matrix
  matrix[n_spp,n_spp] Bmat;
  // diagonal
  Bmat = diag_matrix(Bdiag);
  // off-diagonals
  for(i in 1:n_off) {
    Bmat[rc_off[i,1],rc_off[i,2]] = Boffd[i];
  }
  // expectations
  for(t in 1:n_year) {
    Ey[t] = Zmat * xx[t];
    if (t < 2)
      Ex[t] = Bmat * X0;
    else
      Ex[t] = Bmat * xx[t-1];
  }
  // cov matrix Q
  for(i in 1:n_spp) {
  	QQ[i,i] = SD_proc[id_q[i]];
  }
  for(i in 1:(n_spp-1)) {
  	for(j in (i+1):n_spp) {
  	  QQ[i,j] = 0;
  	  QQ[j,i] = 0;
  	}
  }
  // cov matrix R
  for(i in 1:n_obs) {
  	RR[i,i] = SD_obs[id_r[i]];
  }
  for(i in 1:(n_obs-1)) {
  	for(j in (i+1):n_obs) {
  	  RR[i,j] = 0;
  	  RR[j,i] = 0;
  	}
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
  // for(i in 1:n_q) {
  // 	 SD_proc[i] ~ cauchy(0, 3);
  // }
  SD_proc ~ normal(0, 1);
  // obs SD
  // for(i in 1:n_r) {
  //   SD_obs[i] ~ cauchy(0, 3);
  // }
  SD_obs ~ normal(0, 1);
  // B matrix
  // diagonal
  Bdiag ~ beta(1.05,1.05);
  // Bdiag ~ normal(0.5,10);
  // off-diagonals
  Boffd ~ normal(0,10);
  // likelihood
  xx ~ multi_normal(Ex, QQ);
  yymiss ~ multi_normal(Ey, RR);
}

