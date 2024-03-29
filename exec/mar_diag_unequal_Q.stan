data {
  // scalars
  int<lower=0> n_year;           // # years
  int<lower=0> n_spp;            // # species
  int<lower=0> n_off;            // # non-zero off-diagonals
  int<lower=1> n_q;              // # unique proc SD
  // vectors
  int id_q[n_spp];                 // IDs for proc SD
  // matrices
  matrix[n_spp,n_year] y;        // data
  int<lower=0> rc_off[n_off,2];  // indices of off-diagonals
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
  // vector[n_spp] X0;                    // initial states
  real ymiss[n_na];
}
transformed parameters {
  // cov matrix
  cov_matrix[n_spp] QQ;
  // B matrix
  matrix[n_spp,n_spp] Bmat;
  matrix[n_spp,n_year] yymiss;        // data
  // diagonal
  Bmat = diag_matrix(Bdiag);
  // off-diagonals
  for(i in 1:n_off) {
    Bmat[rc_off[i,1],rc_off[i,2]] = Boffd[i];
  }
  // cov matrix
  for(i in 1:n_spp) {
  	QQ[i,i] = SD_proc[id_q[i]];
  }
  for(i in 1:(n_spp-1)) {
  	for(j in (i+1):n_spp) {
  		QQ[i,j] = 0;
  		QQ[j,i] = 0;
  	}
  }
  // Deal with missing and non-missing values separately
  for(i in 1:n_pos) {
    yymiss[row_indx_pos[i], col_indx_pos[i]] = y[row_indx_pos[i], col_indx_pos[i]];
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
  // X0 ~ normal(0,1);
  // process SD's
  for(i in 1:n_q) {
  	SD_proc[i] ~ cauchy(0, 5);
  }
  // B matrix
  // diagonal
  Bdiag ~ beta(1.05,1.05);
  // off-diagonals
  Boffd ~ normal(0,10);
  // likelihood
  for(t in 2:n_year) {
    col(yymiss,t) ~ multi_normal(Bmat * col(yymiss,t-1), QQ);
  }
}
