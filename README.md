### This is a repository containing code and simulation results for Scheuerell et al. (in review)

The code included is generic, including models for multivariate autogregressive models. We denote these generally as MAR(1) or VAR(1) models when the approach does not include observation or measurement errors, versus MARSS(1) or VARSS(1) models. These models are all stationary, and can be differentiated from models with time-varying parameters included in our other repository, for TVVARSS(1) models https://github.com/nwfsc-timeseries/tvvarsshttps://github.com/nwfsc-timeseries/tvvarss.

### Simulations

The majority of the Bayesian models in our simulations are fit based on a dataframe constructed with the file 'R/make_grid.r'. The Stan files used in our analysis are in the 'exec/' folder, and the code to run these models across all iterations in the grid are in 'R/big_board_sim_fit.r' file.

### Model output and results

Raw results are included in the 'results' folder. Due to file size limits, results have been broken down into a series of files (posterior_summaries_1.rds, posterior_summaries_2.rds, ..., posterior_summaries_9.rds)

Plotting scripts for replicating plots in the paper are in the 'plotting scripts/' folder
Plots are included in the 'plots/' folder

