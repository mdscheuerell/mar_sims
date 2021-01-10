library(dplyr)
library(ggplot2)
library(viridis)
library(ggforce)
grid = readRDS("marss_pars.rds")
grid$pro_sd_label = paste0("pro_sd=",grid$pro_sd)
grid$obs_sd_label = paste0("obs_sd=",grid$obs_sd)
grid$pro_CV_label = paste0("pro_CV=",grid$pro_CV)
grid$obs_CV_label = paste0("obs_CV=",grid$obs_CV)

grid$sd_obs_est = sqrt(grid$R)
grid$sd_pro_est = sqrt(grid$Q)


# make basic plots of
pdf("plots/Figure_S3_obs_error_marss.pdf")
g1 = ggplot(dplyr::filter(grid, b_CV==1, sd_obs_est<1),
  aes(y=sd_obs_est)) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  ylab("Estimated obs error SD") +
  geom_hline(aes(yintercept = obs_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label~ pro_sd_label, scale="free_y")
g1
dev.off()

jpeg("plots/Figure_S3_obs_error_marss.jpeg")
g1
dev.off()

pdf("plots/Figure_S4_pro_error_marss.pdf")
g2 = ggplot(dplyr::filter(grid, b_CV==1,sd_pro_est<1),
  aes(y=sd_pro_est)) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  ylab("Estimated pro error SD") +
  geom_hline(aes(yintercept = pro_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label~ pro_sd_label, scale="free")
g2
dev.off()

jpeg("plots/Figure_S4_pro_error_marss.jpeg")
g2
dev.off()