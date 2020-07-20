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
  int<lower=0> n_pos; // number of non-missing observations
  int<lower=0> row_indx_pos[n_pos];
  int<lower=0> col_indx_pos[n_pos];
  int<lower=0> n_na; // number of missing observations
  int<lower=0> row_indx_na[n_na];
  int<lower=0> col_indx_na[n_na];
  real yy[n_pos];       // data
  real pro_mu; // mean of process variance prior
  real pro_cv; // cv of process variance prior
  real obs_mu; // mean of obs variance prior
  real obs_cv; // cv of obs variance prior
  real b_mu[6]; // prior on mean of offdiagonal
  real b_sd[6]; // prior on sd offdiagonal 
  real b_mu_diag[4];// prior on mean of diagonal
  real b_sd_diag[4];// prior on sd offdiagonal 
}
parameters {
  real<lower=0> SD_obs;
  vector<lower=-1,upper=1>[n_off] Boffd;  // off-diags of B
  vector<lower=0,upper=1>[n_spp] Bdiag;   // diag of B
  real<lower=0> SD_proc;           // proc SD
  matrix[n_spp,n_year] xx;       // states
  vector[n_na] ymiss;  // missing values
}
transformed parameters {
  // obs SD
  // B matrix
  matrix[n_spp,n_spp] Bmat;
  matrix[n_obs,n_year] yymiss;
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
}
model {
  // PRIORS
  // initial state
  xx[,1] ~ normal(0,5);
  // process SD's
  SD_proc ~ normal(pro_mu,pro_mu*pro_cv);
  // obs SD
  SD_obs ~ normal(obs_mu,obs_mu*obs_cv);
  // B diagonal
  //Bdiag ~ normal(0.5,1);
  for(i in 1:4) Bdiag[i] ~ normal(b_mu_diag[i],b_sd_diag[i]);
  // B off-diagonals
  //Boffd ~ normal(0,1);
  for(i in 1:6) Boffd[i] ~ normal(b_mu[i],b_sd[i]);
  // missing obs
  ymiss ~ normal(0,5);
  // LIKELIHOOD
  for(t in 2:n_year) {
    xx[,t] ~ normal(Bmat * xx[,t-1], SD_proc);
    // col(xx,t) ~ multi_normal(Bmat * col(xx,t-1), QQ);
  }
  to_vector(yymiss) ~ normal(to_vector(xx), SD_obs);
}
