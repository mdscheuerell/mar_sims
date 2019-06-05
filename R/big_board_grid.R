
# Proposed simulation framework / datasets

set.seed(123)
replicates = 1 # number of replicates/combination of different parameters

grid = expand.grid("obs_sd" = c(0, 0.2, 0.4, 0.8),
                   "pro_sd" = c(0.1, 0.2, 0.4),
                   "iter" = seq(1,replicates),
                   "frac_missing" = c(0, 0.2, 0.4),
                   "food_web" = c("linear","2_prey","2_basal"),
                   "ts_length" = c(15,30,45))
grid$seed = .Random.seed[grid$iter]

# EW: I propose creating the above grid, but *not necessarily running all of those models*
# It'd help to envision which of the above forms a base case, and we can run models to evaluate 
# each comparison. Maybe set the flag for 'base_case' for the models we want to run = 1
#grid$base_case = 0
# maybe use intermediate values of the above as base case? 
#grid$base_case[which(grid$obs_sd==0.2 & grid$pro_sd==0.2 & grid$frac_missing==0.2 & grid$ts_length==30)] = 1

# Then we can create another flag for whether to run the models, e.g. . This simplifies things a ton, because
# it means we only have to run ~ 500 models instead of 5000 (with 20 replicates). This could be reduced even further 
# if we concentrate just on one food web for these comparisons
#grid$run_model = 0
#grid$run_model[which(grid$base_case==1)] = 1
#grid$run_model[which(grid$pro_sd==0.2 & grid$frac_missing==0.2 & grid$ts_length==30)] = 1 # for obs comparison
#grid$run_model[which(grid$obs_sd==0.2 & grid$frac_missing==0.2 & grid$ts_length==30)] = 1 # for pro comparison
#grid$run_model[which(grid$obs_sd==0.2 & grid$pro_sd==0.2 & grid$ts_length==30)] = 1 # for frac missing comparison
#grid$run_model[which(grid$obs_sd==0.2 & grid$pro_sd==0.2 & grid$frac_missing==0.2)] = 1 # for ts_length comparison

saveRDS(grid, file="grid.rds")