
## Proposed simulation framework / datasets
set.seed(123)
replicates = 100 # number of replicates/combination of different parameters

grid <- expand.grid("obs_sd" = c(0.2, 0.4, 0.8),
  "pro_sd" = c(0.1, 0.2, 0.4),
  "iter" = seq(1,replicates),
  "frac_missing" = c(0),
  "food_web" = c("linear"),
  "ts_length" = c(30),
  "obs_CV"=c(0.1,0.5,1),
  "pro_CV"=c(0.1,0.5,1),
  "b_CV" = c(0.1,0.5,1))

## We don't need all combinations of the above. We need:
## 1. for variance effects, all variance combos with b_CV == 1
grid$keep = 0
grid$keep[which(grid$b_CV==1)] = 1
## 2. for interatction effects, high CV variances (obs_CV, pro_CV==1) and 
grid$keep[which(grid$obs_CV==1 & grid$pro_CV==1)] = 1


grid <- grid[which(grid$keep==1),]
grid$seed <- seq(1,nrow(grid))
saveRDS(grid, file = "grid.rds")

