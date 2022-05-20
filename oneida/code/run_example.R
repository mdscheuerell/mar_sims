
Year = 1977
# prep data
source("oneida/code/04_clean_oneida.R")

#if(length(table(dat_station$station))==3) {
# fit models 
source("oneida/code/05_oneida-fit-1site.R")
source("oneida/code/05_oneida-fit-3site.R")
  
# process
source("oneida/code/06_plot_oneida_results.R")


