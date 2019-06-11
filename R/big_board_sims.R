
##------------------------------
## this is just testing how increasing the magnitude of off-diagonal elements
## impacts our estimates or not
## X[t,1] = B[1,1]*X[t-1,1] + B[1,2]*X[t-1,2] + B[1,3]*X[t-1,3] + B[1,4]*X[t-1,4]
##------------------------------

grid = readRDS("grid.rds")
devtools::install_github("nwfsc-timeseries/tvvarss")
library("tvvarss")
library("rstan")
library("broom")
library("dplyr")

rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

library("here")
stan_dir <- here("exec")

## interaction types
int_types <- c("dd", "td", "bu", "cf")

## number of species
n_species <- 4

for(ii in seq(1,nrow(grid))) {

# set seed
set.seed(grid$seed[ii])

## number of time points
n_time <- grid$ts_length[ii] + 20
## number of initial samples to discard
n_toss <- 20

proc_var_true = grid$pro_sd[ii]^2
obs_var_true <- grid$obs_sd[ii]^2

if(grid$food_web[ii] == "linear") {
  # topo matrix for linear food chain
  B0 <- matrix(list(0),n_species,n_species)
  diag(B0) <- "dd"
  for(i in 1:(n_species-1)) {
    B0[i,i+1] <- "td"
    B0[i+1,i] <- "bu"
  }
  ## row/col indices for off-diag interactions
  rc_off <- do.call(rbind, sapply(int_types[-1], function(x) which(B0==x, arr.ind=T)))
  ## number of off-diag interactions
  n_off <- nrow(rc_off)
  ## true B matrix
  B0_init = diag(0.7,4)
  for(i in 1:3) {
    B0_init[i,i+1]=-0.03
    B0_init[i+1,i] = 0.05
  }
}

## simulate process. var_QX is process error on states.
## var_QB is process var on B -- ignored for static B models.
sim <- simTVVAR(Bt = B0_init,
  topo = B0,
  TT = n_time,
  var_QX = proc_var_true,
  cov_QX = 0,
  var_QB = 0,
  cov_QB = 0)

## true states
xx <- sim$states[,-seq(n_toss+1)]

## observations
yy <- xx + matrix(rnorm(n_species*(n_time-n_toss),0,obs_var_true), n_species, n_time-n_toss)

## number of proc SD's
n_q <- length(unique(proc_var_true))
id_q <- c(1,1,1,1)

## number of obs SD's
n_r <- length(unique(obs_var_true))
id_r <- c(1,1,1,1)

# remove missing data if that's an issue
# algorithm = missing at random
data_to_na = sample(seq(ncol(yy)), size = ncol(yy)*grid$frac_missing[ii])
if(length(data_to_na) > 0) {
  yy[,data_to_na] = NA
}

# indices of positive values - Stan can't handle NAs
row_indx_pos <- matrix(rep(seq_len(nrow(xx)), ncol(xx)), nrow(xx), ncol(xx))[!is.na(xx)]
col_indx_pos <- matrix(sort(rep(seq_len(ncol(xx)), nrow(xx))), nrow(xx), ncol(xx))[!is.na(xx)]
n_pos <- length(row_indx_pos)

row_indx_na <- matrix(rep(seq_len(nrow(xx)), ncol(xx)), nrow(xx), ncol(xx))[is.na(xx)]
col_indx_na <- matrix(sort(rep(seq_len(ncol(xx)), nrow(xx))), nrow(xx), ncol(xx))[is.na(xx)]
n_na <- length(row_indx_na)

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
  rc_off = rc_off
)

## fit model
fit <- stan(file = file.path(stan_dir, "marss_diag_unequal_Q_diag_equal_R.stan"),
            data = dat,
            pars = c("Bmat", "SD_proc", "SD_obs", "xx"),
            control = list(max_treedepth = 25, adapt_delta = 0.99),
            iter = 1000, chains = 3)

# save results to a file
saveRDS(fit, file=paste0("results/raw_models/run_",ii,".rds"))

# also make tidy summary into table -- get 95%, 90% and quartiles
pars = broom::tidy(fit, conf.int=TRUE, conf.level=c(0.95),
  rhat=TRUE, estimate.method = "median")
pars = dplyr::rename(pars, conf_0.025 = conf.low, conf_0.975 = conf.high)
pars2 = broom::tidy(fit, conf.int=TRUE, conf.level=c(0.9),
  rhat=TRUE, estimate.method = "median")
pars$conf_0.05 = pars2$conf.low
pars$conf_0.95 = pars2$conf.high
pars2 = broom::tidy(fit, conf.int=TRUE, conf.level=c(0.5),
  rhat=TRUE, estimate.method = "median")
pars$conf_0.25 = pars2$conf.low
pars$conf_0.75 = pars2$conf.high

pars$iter = ii

if(ii == 1) {
  post_estimates = pars
} else {
  post_estimates = rbind(post_estimates, pars)
}

saveRDS(post_estimates,file="results/posterior_summaries.rds")
}

