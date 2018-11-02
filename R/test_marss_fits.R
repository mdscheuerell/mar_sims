

##------------------------------
## non-diag B; diag & unequal Q
##------------------------------

library("tvvarss")

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library("here")
stan_dir <- here("exec")


## number of species
n_spp <- 4
## number of time points
n_time <- 70
## number of initial samples to discard
n_toss <- 20

## true proc var
# proc_var_true <- rev((seq(n_spp)/10)^2)
proc_var_true <- rev((seq(n_spp)/10))
## true obs var
obs_var_true <- rep(c(2,1)/20,ea=2)

## interaction types
int_types <- c("dd", "td", "bu", "cf")

## topo matrix for linear food chain
B0_lfc <- matrix(list(0),n_spp,n_spp)
diag(B0_lfc) <- "dd"
for(i in 1:(n_spp-1)) {
  B0_lfc[i,i+1] <- "td"
  B0_lfc[i+1,i] <- "bu"
}

## row/col indices for off-diag interactions
rc_off <- do.call(rbind, sapply(int_types[-1], function(x) which(B0_lfc==x, arr.ind=T)))
## number of off-diag interactions
n_off <- nrow(rc_off)

## true B matrix
B0_init <- diag(seq(5,8)/10)
B0_init[1,2] <- runif(1, -0.3, -0.1)
B0_init[2,1] <- runif(1, 0.1, 0.3)
B0_init[2,3] <- runif(1, -0.3, -0.1)
B0_init[3,2] <- runif(1, 0.1, 0.3)
B0_init[3,4] <- runif(1, -0.3, -0.1)
B0_init[4,3] <- runif(1, 0.1, 0.3)
B0_init <- round(B0_init, 2)

## simulate process. var_QX is process error on states.
## var_QB is process var on B -- ignored for static B models.
lfc <- simTVVAR(Bt = B0_init,
  topo = B0_lfc,
  TT = n_time,
  var_QX = proc_var_true,
  cov_QX = 0,
  var_QB = 0,
  cov_QB = 0)

## true states
xx <- lfc$states[,-seq(n_toss+1)]

## observations
yy <- xx + matrix(rnorm(n_spp*(n_time-n_toss),0,obs_var_true), n_spp, n_time-n_toss)

## number of proc SD's
n_q <- length(unique(proc_var_true))
id_q <- seq(n_spp)

## number of obs SD's
n_r <- length(unique(obs_var_true))
id_r <- c(1,1,2,2)
# id_r <- rep(1,4)

## data list for Stan
dat <- list(
  n_year = n_time-n_toss,
  n_species = n_spp,
  n_q = n_q,
  id_q = id_q,
  n_obs = n_spp,
  n_r = n_r,
  id_r = id_r,
  Zmat = diag(n_spp),
  yy = yy,
  n_off = n_off,
  rc_off = rc_off
)

## function to generate inits
initf <- function(chain_id = 1) {
  list(Bdiag = runif(n_spp, 0.5, 0.7), Boffd = runif(n_off, -0.2, 0.2),
       SD_proc=runif(n_q, 0.1,0.2), SD_obs = runif(n_r, 0,0.2))
} 

## list of lists for initial values
n_chains <- 4
inits <- lapply(1:n_chains, function(id) initf(chain_id = id))

## fit model
fit <- stan(file = file.path(stan_dir, "marss_diag_unequal_Q_diag_unequal_R.stan"),
            data = dat,
            pars = c("Bmat", "SD_proc", "SD_obs"),
            init = inits,
            control = list(max_treedepth = 12, adapt_delta = 0.8),
            iter = 10000, thin = 20, chains = n_chains)
fit
## estimated B
round(matrix(summary(fit)$summary[1:(n_spp)^2,"mean"],n_spp,n_spp,byrow=TRUE),2)
## true B
B0_init

traceplot(fit, c("SD_proc","SD_obs"))

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
max_treedepth_by_chain

##----------------------
## unconstrained models
##----------------------

## data list for Stan
dat <- list(
  n_year = n_time-n_toss,
  n_species = n_spp,
  n_q = n_q,
  id_q = id_q,
  n_obs = n_spp,
  n_r = n_r,
  id_r = id_r,
  Zmat = diag(n_spp),
  yy = t(yy),
  n_off = n_off,
  rc_off = rc_off
)

# fit <- stan(file = file.path(stan_dir, "marss_diag_unequal_Q_diag_unequal_R.stan"),
fit <- stan(file = file.path(stan_dir, "marss_uncon_Q_uncon_R.stan"),
            data = dat,
            pars = c("Bmat", "SD_proc", "SD_obs"),
            init = inits,
            control = list(max_treedepth = 15, adapt_delta = 0.9),
            iter = 4000, chains = 2)
fit
## estimated B
round(matrix(summary(fit)$summary[1:(n_spp)^2,"mean"],n_spp,n_spp,byrow=TRUE),2)
## true B
B0_init
