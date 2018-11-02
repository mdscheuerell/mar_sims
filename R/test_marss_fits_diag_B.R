

##------------------------------
## diag B; diag & unequal Q & R
##------------------------------

library("tvvarss")

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library("here")
stan_dir <- here("exec")


## number of species
n_spp <- n_obs <- 4
## number of time points
n_time <- 120
## number of initial samples to discard
n_toss <- 20

## true proc var
# proc_var_true <- rev((seq(n_spp)/10)^2)
proc_var_true <- 2*rev((seq(n_spp)/10))
## true obs var
obs_var_true <- rep(c(2,1)/10,ea=2)
obs_var_true <- rep(1,n_obs)/10

## interaction types
int_types <- c("dd", "td", "bu", "cf")

## topo matrix for linear food chain
B0_lfc <- matrix(list(0),n_spp,n_spp)
diag(B0_lfc) <- "dd"

## true B matrix
B0_init <- diag(seq(5,8)/10)
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
# id_r <- c(1,1,2,2)
id_r <- rep(1,4)

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
  yy = yy
)

## function to generate inits
initf <- function(chain_id = 1) {
  list(Bdiag = runif(n_spp, 0.5, 0.7),
       SD_proc=runif(n_q, 0.1,0.2), SD_obs = runif(n_r, 0.1 ,0.2))
} 

## list of lists for initial values
n_chains <- 4
inits <- lapply(1:n_chains, function(id) initf(chain_id = id))

## fit model
fit <- stan(file = file.path(stan_dir, "marss_diag_B_diag_unequal_Q_diag_unequal_R.stan"),
            data = dat,
            pars = c("Bdiag", "SD_proc", "SD_obs"),
            init = inits,
            control = list(max_treedepth = 15, adapt_delta = 0.95),
            iter = 6000, warmup = 4000, thin = 8, chains = n_chains)
fit
## estimated B
round(diag(summary(fit)$summary[1:(n_spp),"mean"]),2)
## true B
B0_init

traceplot(fit, c("SD_proc","SD_obs"))

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
max_treedepth_by_chain

