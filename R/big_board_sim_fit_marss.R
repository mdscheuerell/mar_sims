
##-------------------
## initialize & load
if(!require("tvvarss")) {
  install.packages("tvvarss")
  library("tvvarss")
}
if(!require("MARSS")) {
  install.packages("MARSS")
  library("MARSS")
}
if(!require("here")) {
  install.packages("here")
  library("here")
}

## set dir locations
stan_dir <- here("exec")
res_dir <- here("results")
raw_dir <- here("results/raw_models")

## load grid of sim options
grid <- readRDS("grid.rds")

##-----------------
## specify options
##-----------------

## interaction types
int_types <- c("dd", "td", "bu", "cf")

## number of initial samples to discard
n_toss <- 20

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

##-----------
## sim & fit
##-----------

for(ii in 1:nrow(grid)) {

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

  fit = MARSS(y = yy, model = list("R"="diagonal and equal","Q"="diagonal and equal", "U"="zero",
                                   "B"=matrix(list("b11","b21",0,0,"b12","b22","b32",0,0,"b23","b33","b43",0,0,"b34","b44"),4,4)),
              control=list(maxit=1000,conv.test.slope.tol=1000),silent=TRUE,fun.kf="MARSSkfas")
  
  ## save raw results to a file
  saveRDS(fit, file = file.path(raw_dir, paste0("marss_", ii, ".rds")))

}
grid$id = seq(1,nrow(grid))
marss_pars = matrix(nrow=nrow(grid), ncol = 13)
for(ii in 1:nrow(grid)) {
  fit = readRDS(file.path(raw_dir, paste0("marss_", ii, ".rds")))
  pars = c(fit$par$B, fit$par$Q, fit$par$Q)
  marss_pars[ii,1] = ii
  marss_pars[ii,2:13] = c(fit$par$B, fit$par$R, fit$par$Q)
}
colnames(marss_pars) = c("id","b11","b21","b12","b22","b32","b23","b33","b43","b34","b44","R","Q")
marss_pars = cbind(grid, marss_pars)
saveRDS(marss_pars,"marss_pars.rds")



#g = group_by(posterior_summaries_parII[grep("Bmat",rownames(posterior_summaries_parII)),], iter) %>%
#  summarize(m = max(Rhat,na.rm=T))

#g1 = group_by(posterior_summaries_parII[grep("SD_proc",rownames(posterior_summaries_parII)),], iter) %>%
#  summarize(m = max(Rhat,na.rm=T))

#g2 = group_by(posterior_summaries_parII[grep("SD_obs",rownames(posterior_summaries_parII)),], iter) %>%
#  summarize(m = max(Rhat,na.rm=T))

#g3 = group_by(posterior_summaries_parII[grep("xx",rownames(posterior_summaries_parII)),], iter) %>%
#  summarize(m = max(Rhat,na.rm=T))
