
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
library(ggmcmc)
library(dplyr)
## set stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## set dir locations
stan_dir <- here("oneida/models")
res_dir <- here("oneida/output")
raw_dir <- here("oneida/output/raw_models")

## interaction types
int_types <- c("dd", "td", "bu", "cf")

## stan options
stan_model <- "marss_diag_unequal_Q_diag_equal_prior_R.stan"#"marss_diag_unequal_Q_diag_equal_prior_R_2survey.stan"
stan_ctrl <- list(max_treedepth = 15, adapt_delta = 0.99)
stan_mcmc <- list(iter = 3000, warmup = 2000, chains = 3, thin = 1, refresh = 0)
## interaction matrices (B)
## index 1 is bottom of food web; 4 is top

dat = readRDS("oneida/data/cleaned_data_agg.rds")
dat = as.matrix(dat[,-1])
# The phytoplankton community is largely comprised of the five
# major taxonomic groups typical of North Temperate Zone lakes including Bacillariophyta,
# Cryptophyta, Chrysophyta, Chlorophyta, and Cyanobacteria. In most years, diatoms are the
# dominant taxon in early spring and fall, small flagellated cryptophytes and chlorophytes
# dominate the assemblage (though at low densities) during the clear water phase, with
# cyanobacterial blooms taking place between July and October, including in 2015

## initialize table of posterior summaries
post_estimates <- NULL

## number of species
n_species <- 4

## observations
  yy <- dat
  # yy <- t(scale(t(yy)))
  #yy <- t(scale(t(yy)))
  n_year <- ncol(yy)
  ## number of proc SD's
  n_q <- 1
  id_q <- c(1,1,1,1)

  ## number of obs SD's
  n_r <- 1
  id_r <- c(1,1,1,1)

  B_topo <- list("dd", 0, "td", "td", # BlueGreen
                 0, "dd", "td", "td", # green 
                 "bu", "bu", "dd", "bu", # Daphnia
                 0,    0, 0, "dd") # Diatoms
  topo <- matrix(B_topo, n_species, n_species, byrow = TRUE)
  
  ## row/col indices for off-diagonals
  rc_off <- do.call(rbind, sapply(int_types[-1], function(x) which(topo == x, arr.ind = TRUE)))
  ## number of non-zero off-diagonals
  n_off <- nrow(rc_off)
  
  ## indices for non-NA values
  row_indx_pos <- matrix(rep(seq_len(nrow(yy)), ncol(yy)), nrow(yy), ncol(yy))[!is.na(yy)]
  col_indx_pos <- matrix(sort(rep(seq_len(ncol(yy)), nrow(yy))), nrow(yy), ncol(yy))[!is.na(yy)]
  n_pos <- length(row_indx_pos)

  ## data without NA
  yy <- yy[!is.na(yy)]

  ## data list for Stan
  dat <- list(
    n_year = n_year,
    n_spp = n_species,
    n_q = n_q,
    id_q = id_q,
    n_obs = n_species,
    n_r = n_r,
    id_r = id_r,
    Zmat = diag(n_species),
    yy = yy,
    #yy2 = yy2,
    n_off = n_off,
    rc_off = rc_off,
    n_pos= n_pos,
    row_indx_pos = row_indx_pos,
    col_indx_pos = col_indx_pos,
    # n_na = n_na,
    # row_indx_na = row_indx_na,
    # col_indx_na = col_indx_na,
    pro_mu = 0.2,
    pro_cv = 1,
    obs_mu = 0.1,
    obs_cv = 1,
    b_mu = c(-0.1,-0.1,-0.1,-0.1,0.1,0.1,0),
    b_sd = rep(0.03, 7),
    b_mu_diag = rep(0.5, 4),
    b_sd_diag = rep(0.3, 4)
  )

  ## fit model
  stan_model <- "cbfs_1site.stan"
  fit <- stan(file = file.path(stan_dir, stan_model),
                  data = dat,
                  pars = c("Bmat", "SD_proc", "SD_obs","xx"),
                  control = stan_ctrl,
                  iter = stan_mcmc$iter,
                  warmup = stan_mcmc$warmup,
                  chains = stan_mcmc$chains,
                  thin = stan_mcmc$thin)#,
                  #refresh = stan_mcmc$refresh)

  saveRDS(fit, "oneida/output/fitted_model_1site.rds")

