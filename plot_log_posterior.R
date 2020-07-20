
# SK - change this next line
this_batch = 1 # 1, 2, 3

n_batch = 3 # max - 3 by default

##-------------------
## initialize & load
##-------------------

if(!require("tvvarss")) {
  devtools::install_github("nwfsc-timeseries/tvvarss")
  library("tvvarss")
}
if(!require("rstan")) {
  install.packages("rstan")
  library("rstan")
}
if(!require("here")) {
  install.packages("here")
  library("here")
}

## set stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## set dir locations
stan_dir <- here("exec")
res_dir <- here("results")
raw_dir <- here("results/raw_models")

## load grid of sim options
grid <- readRDS("grid.rds")
grid$batch = sort(rep(1:n_batch, ceiling(nrow(grid)/n_batch)))[1:nrow(grid)]

##-----------------
## specify options
##-----------------

## interaction types
int_types <- c("dd", "td", "bu", "cf")

## number of initial samples to discard
n_toss <- 20

## stan options
stan_model <- "marss_diag_unequal_Q_diag_equal_prior_R.stan"
stan_ctrl <- list(max_treedepth = 25, adapt_delta = 0.999)
stan_mcmc <- list(iter = 3000, warmup = 2000, chains = 3, thin = 1, refresh = 0)

## interaction matrices (B)
## index 1 is bottom of food web; 4 is top

## linear food chain (eg, phyto -> zoop -> minnows -> bass)
topo_lfc <- list("dd", "td",    0,    0,
                 "bu", "dd", "td",    0,
                 0, "bu", "dd", "td",
                 0,    0, "bu", "dd")
B0_lfc <- c(0.5, -0.1,  0.0,  0.0,
            0.3,  0.6, -0.2,  0.0,
            0.0,  0.2,  0.7, -0.3,
            0.0,  0.0,  0.1,  0.8)
B0_mat = matrix(B0_lfc,4,4, byrow=TRUE)

##-------
## setup
##-------

## initialize table of posterior summaries
post_estimates <- NULL

## number of species
n_species <- sqrt(length(B0_lfc))

## initial values
init_vals <- function(chain_id = 1, n_off, n_species, n_time, n_na) {
  list(SD_obs = runif(1),
       SD_proc = runif(1),
       Boffd = runif(n_off, -0.5, 0.5),
       Bdiag = runif(n_species),
       xx = matrix(rnorm(n_species*n_time), n_species, n_time),
       # init_state = rnorm(4, 0, 5),
       # Z_proc = matrix(rnorm(n_species*n_time), n_species, n_time),
       ymiss = rnorm(n_na)
  )
}
# init_vals <- function(chain_id = 1, n_off, n_species, n_time, n_na) {
#   list(var_SD = runif(1),
#        cov_SD = runif(1,-1,0),
#        mean_SD = runif(1),
#        SD_vec = runif(2),
#        Boffd = runif(n_off, -0.5, 0.5),
#        Bdiag = runif(n_species),
#        xx = matrix(rnorm(n_species*n_time), n_species, n_time),
#        ymiss = rnorm(n_na))
# }


##-----------
## sim & fit
##-----------

ii = which(grid$obs_CV==1 & grid$pro_CV==1 & grid$b_CV==1 & grid$obs_sd==0.2 & grid$pro_sd==0.2)[1]

  ## set seed
  set.seed(grid$seed[ii])
  
  ## number of time points
  n_time <- grid$ts_length[ii] + n_toss
  
  ## proc & obs var
  proc_sd_true <- grid$pro_sd[ii]
  obs_sd_true <- grid$obs_sd[ii]
  
  ## get correct B matrix & topology
  if(grid$food_web[ii] == "linear") {
    B_vals <- B0_lfc
    B_topo <- topo_lfc
  } else if(grid$food_web[ii] == "2_prey") {
    B_vals <- B0_int
    B_topo <- topo_int
  } else {
    B_vals <- B0_bas
    B_topo <- topo_bas
  }
  
  topo <- matrix(B_topo, n_species, n_species, byrow = TRUE)
  
  ## row/col indices for off-diagonals
  rc_off <- do.call(rbind, sapply(int_types[-1], function(x) which(topo == x, arr.ind = TRUE)))
  ## number of non-zero off-diagonals
  n_off <- nrow(rc_off)
  
  ## simulate process. var_QX is process error on states.
  ## var_QB is process var on B -- ignored for static B models.
  sim <- simTVVAR(Bt = B0_mat,
                  topo = topo,
                  TT = n_time,
                  var_QX = proc_sd_true^2,
                  cov_QX = 0,
                  var_QB = 0,
                  cov_QB = 0)
  
  ## true states
  xx <- sim$states[,-seq(n_toss+1)]
  
  ## observations
  yy <- xx + matrix(rnorm(n_species*(n_time-n_toss), 0, obs_sd_true), n_species, n_time-n_toss)
  # yy <- t(scale(t(yy)))
  
  ## number of proc SD's
  n_q <- length(unique(proc_sd_true))
  id_q <- c(1,1,1,1)
  
  ## number of obs SD's
  n_r <- length(unique(obs_sd_true))
  id_r <- c(1,1,1,1)
  
  ## remove missing data if that's an issue
  ## algorithm = missing at random
  data_to_na <- sample(seq(ncol(yy)), size = ncol(yy)*grid$frac_missing[ii])
  if(length(data_to_na) > 0) {
    yy[,data_to_na] <- NA
  }
  
  ## indices for non-NA values
  row_indx_pos <- matrix(rep(seq_len(nrow(yy)), ncol(yy)), nrow(yy), ncol(yy))[!is.na(yy)]
  col_indx_pos <- matrix(sort(rep(seq_len(ncol(yy)), nrow(yy))), nrow(yy), ncol(yy))[!is.na(yy)]
  n_pos <- length(row_indx_pos)
  
  ## indices for NA values
  row_indx_na <- matrix(rep(seq_len(nrow(yy)), ncol(yy)), nrow(yy), ncol(yy))[is.na(yy)]
  col_indx_na <- matrix(sort(rep(seq_len(ncol(yy)), nrow(yy))), nrow(yy), ncol(yy))[is.na(yy)]
  n_na <- length(row_indx_na)
  
  ## data without NA
  yy <- yy[!is.na(yy)]
  
  ## get mean of B elemenets
  b_mu = rep(0, nrow(rc_off))
  b_sd = rep(0, nrow(rc_off))
  for(jj in 1:nrow(rc_off)) {
    b_mu[jj] = B0_mat[rc_off[jj,1],rc_off[jj,2]]
    b_sd[jj] = grid$b_CV[ii]*b_mu[jj]
  }
  ## data list for Stan
  dat <- list(
    n_year = n_time-n_toss,
    n_spp = n_species,
    n_q = n_q,
    id_q = id_q,
    n_obs = n_species,
    n_r = n_r,
    id_r = id_r,
    Zmat = diag(n_species),
    yy = yy,
    n_off = n_off,
    rc_off = rc_off,
    n_pos= n_pos,
    row_indx_pos = row_indx_pos,
    col_indx_pos = col_indx_pos,
    n_na = n_na,
    row_indx_na = row_indx_na,
    col_indx_na = col_indx_na,
    pro_mu = grid$pro_sd[ii],
    pro_cv = grid$pro_CV[ii],
    obs_mu = grid$obs_sd[ii],
    obs_cv = grid$obs_CV[ii],
    b_mu = b_mu,
    b_sd = abs(b_sd),
    b_mu_diag = diag(B0_mat),
    b_sd_diag = diag(B0_mat)*grid$b_CV[ii]
  )
  
  ## initial values
  init_ll <- lapply(1:stan_mcmc$chains,
                    function(id) init_vals(chain_id = id,
                                           n_off = n_off,
                                           n_species = n_species,
                                           n_time = grid$ts_length[ii],
                                           n_na = n_na)
  )

# loop through ~ 100 random seeds, using the same dataset above
set.seed(123)

  ## fit model, without passing init values in
  fit <- try(stan(file = file.path(stan_dir, stan_model),
                  data = dat,
                  pars = c("Bmat", "SD_proc", "SD_obs"),
                  control = stan_ctrl,
                  iter = 30000,
                  warmup = 5000,
                  chains = stan_mcmc$chains,
                  thin = stan_mcmc$thin,
                  refresh = stan_mcmc$refresh),
             silent=TRUE)
saveRDS(fit, "fit.rds")

fit = readRDS("fit.rds")

pdf("plots/pairs_fit.pdf")
  pairs(fit)
dev.off()

fit = as.data.frame(fit)
fitdf = fit[,c("Bmat[1,1]","Bmat[1,2]","SD_proc","SD_obs","lp__")]
names(fitdf) = c("B11","B12","sd_proc","sd_obs","lp")

g1 = ggplot(fitdf, aes(x=sd_obs, y=sd_proc, col=lp)) + 
  geom_hex(bins=50) + scale_fill_viridis_c() + 
  xlab(expression(paste("Obs. error ",sigma))) + 
  ylab(expression(paste("Pro. error ",sigma))) + 
  labs(fill = "Frequency")

g2 = ggplot(fitdf, aes(sd_obs, sd_proc, z=lp)) + 
  stat_summary_hex(aes(z = lp), bins = 50, fun = mean) + scale_fill_viridis_c() + 
  xlab(expression(paste("Obs. error ",sigma))) + 
  ylab(expression(paste("Pro. error ",sigma))) + 
  labs(fill = "Log posterior")

pdf("plots/variance_tradeoff.pdf")
gridExtra::grid.arrange(g1,g2,nrow=1)
dev.off()

# subtract off the prior on the variance parameters
fitdf$ll = fitdf$lp - dnorm(fitdf$sd_obs, grid$obs_sd[ii],grid$obs_sd[ii],log=TRUE) - dnorm(fitdf$sd_pro, grid$pro_sd[ii],grid$pro_sd[ii],log=TRUE)
ggplot(fitdf, aes(sd_obs, sd_proc, z=ll)) + 
  stat_summary_hex(aes(z = ll), bins = 30, fun = mean) + scale_fill_viridis_c() + 
  xlab(expression(paste("Obs. error ",sigma))) + 
  ylab(expression(paste("Pro. error ",sigma))) + 
  labs(fill = "Log likelihood")


g3 = ggplot(fitdf, aes(x=B11, y=B12, col=lp)) + 
  geom_hex(bins=50) + scale_fill_viridis_c() + 
  xlab("B[1,1]") + 
  ylab("B[1,2]") + 
  labs(fill = "Frequency")

g4 = ggplot(fitdf, aes(B11, B12, z=lp)) + 
  stat_summary_hex(aes(z = lp), bins = 50, fun = mean) + scale_fill_viridis_c() + 
  xlab("B[1,1]") + 
  ylab("B[1,2]") + 
  labs(fill = "Log posterior")

pdf("plots/B11_B12_tradeoff.pdf")
gridExtra::grid.arrange(g3,g4,nrow=1)
dev.off()
