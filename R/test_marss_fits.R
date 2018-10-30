

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
n_time <- 120
## number of initial samples to discard
n_toss <- 20

## true proc var
proc_var_true <- rev((seq(n_spp)/10)^2)
## true obs var
obs_var_true <- rep(c(1,2)/10,ea=2)^2

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

## number of proc SD's
n_r <- length(unique(obs_var_true))
id_r <- c(1,1,2,2)

## data list for Stan
dat <- list(
  n_time = n_time-n_toss,
  n_spp = n_spp,
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

## fit model
fit <- stan(file = file.path(stan_dir, "marss_diag_unequal_Q_diag_unequal_R.stan"),
            data = dat,
            pars = c("Bmat", "SD_proc", "SD_obs"),
            control = list(max_treedepth = 15),
            iter = 4000, chains = 2)
fit
## estimated B
round(matrix(summary(fit)$summary[1:(n_spp)^2,"mean"],n_spp,n_spp,byrow=TRUE),2)
## true B
B0_init
