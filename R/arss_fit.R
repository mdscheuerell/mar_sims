
library("rstan")
library("here")

## set stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## set dir locations
stan_dir <- here("exec")

## stan options
stan_model <- "arss_nc.stan"
stan_ctrl <- list(max_treedepth = 15, adapt_delta = 0.999)
stan_mcmc <- list(iter = 6000, warmup = 5000, chains = 4, thin = 10, refresh = 0)

## initial values
init_vals <- function(chain_id = 1, n_time) {
  list(SD_obs = runif(1, 0.1, 1),
       SD_proc = runif(1, 0.1, 1),
       phi = runif(1, 0.1, 0.9),
       init_state = rnorm(1),
       zz = rnorm(n_time)
  )
}

## AR(1) params for simulation
nn <- 50
phi <- 0.7
SD_proc <- sqrt(0.4)
SD_obs <- sqrt(0.1)
init_state <- 0

## init states & proc errors
xx <- ww <- rnorm(nn, 0, SD_proc)

## states
xx[1] <- phi * init_state + ww[1]
for(t in 2:nn) {
  xx[t] <- phi * xx[t-1] + ww[t]
}

## obs data
yy <- xx + rnorm(nn, 0, SD_obs)


## data list for Stan
dat <- list(n_year = nn, yy = yy)

## initial values
init_ll <- lapply(1:stan_mcmc$chains,
                  function(id) init_vals(chain_id = id, n_time = nn)
                  )

## fit the model
fit <- stan(file = file.path(stan_dir, stan_model),
            data = dat,
            pars = c("phi", "SD_proc", "SD_obs"),
            control = stan_ctrl,
            init = init_ll,
            iter = stan_mcmc$iter,
            warmup = stan_mcmc$warmup,
            chains = stan_mcmc$chains,
            thin = stan_mcmc$thin,
            refresh = stan_mcmc$refresh)


## SHIELD YOUR EYES BEFORE EXECUTING THESE COMMANDS

## examine fits
print(fit, pars = c("SD_proc", "SD_obs", "phi"))

traceplot(fit, pars = c("SD_proc", "SD_obs", "phi"), inc_warmup = FALSE, nrow = 3)

shinystan::launch_shinystan(fit)
