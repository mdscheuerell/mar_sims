## load MARSS
library(MARSS)

## set directory
dat_dir <- here::here("results")

## load Eric's saved states & obs
load(file.path(dat_dir, "Eric_sim_data.Rdata"))

## define some parameters
## number of observations
nn <- 2
## number of states
mm <- 4

## define interaction matrix (B)
BB <- matrix(list(0), mm, mm)
diag(BB) <- paste0("b", seq(4), seq(4))
for(i in 1:3) {
  BB[i,i+1] <- paste0("b", i+1, i)
  BB[i+1,i] <- paste0("b", i, i+1)
}

## define proc cov amtrix (Q)
QQ <- matrix(list(0), mm, mm)
diag(QQ) <- "q"

#---------------------
# model with 1 obs ts
#---------------------

## define obs cov matrix (R)
R1 <- matrix(list(0), mm, mm)
diag(R1) <- "r"

## define state-space mapping
Z1 <- diag(mm)

## define MARSS model for 1 obs
mod_list_1 <- list(
  B = BB,
  U = matrix(0, mm, 1),
  Q = QQ,
  Z = Z1,
  A = matrix(0, mm, 1),
  R = R1
)

## fit MARSS model to each of the 2 obs ts
mod_obs_1a <- MARSS(yy, mod_list_1)
mod_obs_1b <- MARSS(yy2, mod_list_1)

## estimated B
round(coef(mod_obs_1a, type = "matrix")$B, 2)
#      [,1]  [,2]  [,3]  [,4]
# [1,] 0.51 -0.12  0.00  0.00
# [2,] 0.28  0.65 -0.24  0.00
# [3,] 0.00  0.20  0.67 -0.32
# [4,] 0.00  0.00  0.12  0.80

round(coef(mod_obs_1b, type = "matrix")$B, 2)
#      [,1]  [,2]  [,3]  [,4]
# [1,] 0.49 -0.11  0.00  0.00
# [2,] 0.27  0.61 -0.27  0.00
# [3,] 0.00  0.23  0.68 -0.31
# [4,] 0.00  0.00  0.14  0.80

## bootstrapped CI's
MARSSparamCIs(mod_obs_1a)
#         ML.Est Std.Err  low.CI   up.CI
# R.r      0.165  0.0169  0.1323  0.1984
# B.b11    0.515  0.0609  0.3954  0.6343
# B.b12    0.276  0.0583  0.1613  0.3899
# B.b21   -0.122  0.0388 -0.1985 -0.0465
# B.b22    0.653  0.0422  0.5699  0.7352
# B.b23    0.203  0.0400  0.1244  0.2814
# B.b32   -0.238  0.0351 -0.3064 -0.1688
# B.b33    0.672  0.0401  0.5929  0.7503
# B.b34    0.121  0.0363  0.0501  0.1923
# B.b43   -0.317  0.0418 -0.3991 -0.2351
# B.b44    0.798  0.0436  0.7128  0.8839
# Q.q      0.162  0.0207  0.1212  0.2021
# x0.X.Y1 -1.397  1.0065 -3.3701  0.5754
# x0.X.Y2  0.645  0.8420 -1.0057  2.2949
# x0.X.Y3 -0.534  0.7760 -2.0551  0.9867
# x0.X.Y4  0.858  0.6285 -0.3735  2.0902

MARSSparamCIs(mod_obs_1b)
#         ML.Est Std.Err  low.CI   up.CI
# R.r      0.176  0.0179  0.1411  0.2114
# B.b11    0.492  0.0634  0.3679  0.6166
# B.b12    0.267  0.0598  0.1500  0.3845
# B.b21   -0.107  0.0414 -0.1884 -0.0263
# B.b22    0.609  0.0449  0.5213  0.6975
# B.b23    0.234  0.0431  0.1492  0.3180
# B.b32   -0.265  0.0364 -0.3364 -0.1939
# B.b33    0.680  0.0409  0.6003  0.7607
# B.b34    0.136  0.0367  0.0637  0.2076
# B.b43   -0.311  0.0415 -0.3926 -0.2301
# B.b44    0.803  0.0435  0.7178  0.8881
# Q.q      0.162  0.0217  0.1199  0.2048
# x0.X.Y1 -1.390  1.0801 -3.5069  0.7272
# x0.X.Y2  1.001  0.9062 -0.7747  2.7775
# x0.X.Y3 -0.469  0.7639 -1.9666  1.0280
# x0.X.Y4  0.532  0.6315 -0.7053  1.7701


#---------------------
# model with 2 obs ts
#---------------------

## define obs cov matrix (R)
R2 <- matrix(list(0), nn*mm, nn*mm)
diag(R2) <- "r"

## define state-space mapping
Z2 <- rbind(diag(mm), diag(mm))

## define MARSS model for 1 obs
mod_list_2 <- list(
  B = BB,
  U = matrix(0, mm, 1),
  Q = QQ,
  Z = Z2,
  A = matrix(0, mm*nn, 1),
  R = R2
)

## fit MARSS model to both of the obs ts
mod_obs_2 <- MARSS(rbind(yy, yy2), mod_list_2)

## estimated B
round(coef(mod_obs_2, type = "matrix")$B, 2)
#      [,1]  [,2]  [,3]  [,4]
# [1,] 0.47 -0.10  0.00  0.00
# [2,] 0.26  0.62 -0.23  0.00
# [3,] 0.00  0.21  0.66 -0.32
# [4,] 0.00  0.00  0.11  0.79

## bootstrapped CI's
MARSSparamCIs(mod_obs_2)
#       ML.Est Std.Err low.CI   up.CI
# R.r    0.158 0.00490  0.148  0.1671
# B.b11  0.473 0.04744  0.380  0.5661
# B.b12  0.257 0.04794  0.163  0.3512
# B.b21 -0.103 0.03595 -0.174 -0.0329
# B.b22  0.620 0.03662  0.548  0.6917
# B.b23  0.206 0.03669  0.134  0.2776
# B.b32 -0.234 0.03222 -0.298 -0.1713
# B.b33  0.662 0.03407  0.595  0.7288
# B.b34  0.112 0.03262  0.048  0.1758
# B.b43 -0.321 0.03730 -0.394 -0.2478
# B.b44  0.785 0.03589  0.715  0.8554
# Q.q    0.179 0.00939  0.160  0.1969
# x0.X1 -1.327 1.00141 -3.290  0.6358
# x0.X2  0.875 0.81115 -0.715  2.4650
# x0.X3 -0.453 0.72839 -1.881  0.9745
# x0.X4  0.843 0.59889 -0.330  2.0173

