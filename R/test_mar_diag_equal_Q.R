

##------------------------------
## new stuff for non-diagonal B
##------------------------------

library("tvvarss")

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library("here")
stan_dir <- here("exec")

n_species <- 4

n_year <- 60

## true proc var
proc_var_true <- 0.1^2

## interaction types
int_types <- c("dd", "td", "bu", "cf")

## topo matrix for linear food chain
B0_lfc <- matrix(list(0),n_species,n_species)
diag(B0_lfc) <- "dd"
for(i in 1:(n_species-1)) {
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
  TT = n_year,
  var_QX = proc_var_true,
  cov_QX = 0,
  var_QB = 0,
  cov_QB = 0)

## true states
xx <- lfc$states[,-(1:21)]

## data list for Stan
dat <- list(
n_year = n_year-20,
n_spp = n_species,
n_off = n_off,
y = xx,
rc_off
)

## fit model
fit <- stan(file = file.path(stan_dir,"mar_diag_equal_Q.stan"),
            data = dat, pars = c("Bmat", "sigma"),
            iter = 3000, chains = 2)
## estimated B
round(matrix(summary(fit)$summary[1:16,"mean"],n_species,n_species,byrow=TRUE),2)
## true B
B0_init
            