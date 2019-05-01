
##------------------------------
## this is just testing how increasing the magnitude of off-diagonal elements
## impacts our estimates or not
## X[t,1] = B[1,1]*X[t-1,1] + B[1,2]*X[t-1,2] + B[1,3]*X[t-1,3] + B[1,4]*X[t-1,4]
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
proc_var_true <- 0.03
## true obs var
obs_var_true <- 0.1

## interaction types
int_types <- c("dd", "td", "bu", "cf")

## topo matrix for linear food chain
B0_lfc <- matrix(list(0),n_spp,n_spp)
diag(B0_lfc) <- "dd"
B0_lfc[1,2] = "td"

## row/col indices for off-diag interactions
rc_off <- do.call(rbind, sapply(int_types[-1], function(x) which(B0_lfc==x, arr.ind=T)))
## number of off-diag interactions
n_off <- nrow(rc_off)

sims = expand.grid("B_off" = -c(0.025,0.05,0.075,0.1,0.125,0.15), "sims"=1:50,
  "B_diag" = c(0.5,0.75,1))
sims$B11_25 = NA
sims$B11_50 = NA
sims$B11_75 = NA
sims$B12_25 = NA
sims$B12_50 = NA
sims$B12_75 = NA
sims$B22_25 = NA
sims$B22_50 = NA
sims$B22_75 = NA
sims$B33_25 = NA
sims$B33_50 = NA
sims$B33_75 = NA
sims$B44_25 = NA
sims$B44_50 = NA
sims$B44_75 = NA

for(ii in 1:nrow(sims)) {

set.seed(sims$sims[ii])
## true B matrix
B0_init <- diag(sims$B_diag[ii],4)
B0_init[1,2] <- sims$B_off[ii]

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
id_q <- c(1,1,1,1)

## number of obs SD's
n_r <- length(unique(obs_var_true))
id_r <- c(1,1,1,1)
# id_r <- rep(1,4)

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

## fit model
fit <- stan(file = file.path(stan_dir, "marss_diag_unequal_Q_diag_unequal_R.stan"),
            data = dat,
            pars = c("Bmat", "SD_proc", "SD_obs"),
            control = list(max_treedepth = 25, adapt_delta = 0.99),
            iter = 2000, chains = 3)

B = extract(fit,"Bmat",permuted=TRUE)$Bmat
sims$B11_25[ii] = quantile(B[,1,1],0.25)
sims$B11_50[ii] = quantile(B[,1,1],0.5)
sims$B11_75[ii] = quantile(B[,1,1],0.75)
sims$B12_25[ii] = quantile(B[,1,2],0.25)
sims$B12_50[ii] = quantile(B[,1,2],0.5)
sims$B12_75[ii] = quantile(B[,1,2],0.75)
sims$B22_25[ii] = quantile(B[,2,2],0.25)
sims$B22_50[ii] = quantile(B[,2,2],0.5)
sims$B22_75[ii] = quantile(B[,2,2],0.75)
sims$B33_25[ii] = quantile(B[,3,3],0.25)
sims$B33_50[ii] = quantile(B[,3,3],0.5)
sims$B33_75[ii] = quantile(B[,3,3],0.75)
sims$B44_25[ii] = quantile(B[,4,4],0.25)
sims$B44_50[ii] = quantile(B[,4,4],0.5)
sims$B44_75[ii] = quantile(B[,4,4],0.75)

}

save(sims,file="test_off_diagonal_b_matrix.Rdata")

ggplot(sims, aes(B_off, B12_50-B_off,group=B_off)) + geom_boxplot()

ggplot(sims, aes(B_off, B22_50,group=B_off)) + geom_boxplot() +
  facet_wrap(~B_diag)

pdf("test_off_diagonal_b.pdf")
ggplot(sims, aes(B_off, abs(B12_50-B_off)/(B_off),group=B_off)) + geom_boxplot() +
  facet_wrap(~B_diag) + xlab("True B[1,2]") + ylab("Percent abs. error")
dev.off()

# we can also plot individual datasets
ggplot(sims,
  aes(B_off, B12_50, group=sims, col=sims)) +
  geom_line()
