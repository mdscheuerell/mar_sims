library(dplyr)
library(ggplot2)

grid = readRDS("marss_pars.rds")
grid$pro_sd_label = paste0("pro_sd=",grid$pro_sd)
grid$obs_sd_label = paste0("obs_sd=",grid$obs_sd)
grid$pro_CV_label = paste0("pro_CV=",grid$pro_CV)
grid$obs_CV_label = paste0("obs_CV=",grid$obs_CV)

grid$sd_obs_est = sqrt(grid$R)
grid$sd_pro_est = sqrt(grid$Q)
  
# make basic plots of
pdf("plots/Figure_S3.pdf")
g1 = ggplot(dplyr::filter(grid, b_CV==1,sd_obs_est < 1),
       aes(x = as.factor(""), y=sd_obs_est)) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  ylab(expression(paste("Estimated observation error ",sigma))) +
  geom_hline(aes(yintercept = obs_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label~ pro_sd_label) + 
  xlab("Estimate")+
  theme_bw() + 
  theme(legend.position='none',strip.background = element_rect(color="black",fill="white"))
g1
dev.off()

jpeg("plots/Figure_S3.jpeg")
g1
dev.off()

pdf("plots/Figure_S4.pdf")
g2 = ggplot(dplyr::filter(grid, b_CV==1, sd_pro_est<1),
  aes(x = as.factor(""), y=sd_pro_est)) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  ylab(expression(paste("Estimated process error ",sigma))) +
  geom_hline(aes(yintercept = pro_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label~ pro_sd_label)+
  xlab("Estimate")+
  theme_bw() + 
  theme(legend.position='none',strip.background = element_rect(color="black",fill="white"))
g2
dev.off()

jpeg("plots/Figure_S4.jpeg")
g2
dev.off()
