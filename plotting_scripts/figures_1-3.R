library(dplyr)
library(ggplot2)
grid = readRDS("grid.rds")
post = rbind(readRDS("results/posterior_summaries_t_1.rds"),
             readRDS("results/posterior_summaries_t_2.rds"),
             readRDS("results/posterior_summaries_t_3.rds"))
post$par = rownames(post)

sd_obs_locs = grep("SD_obs",rownames(post))
grid$sd_obs_est = post$mean[sd_obs_locs]
sd_pro_locs = grep("SD_pro",rownames(post))
grid$sd_pro_est = post$mean[sd_pro_locs]

# labels
grid$pro_sd_label = factor(grid$pro_sd, 
  labels = c("sigma[pro] == 0.1","sigma[pro] == 0.2","sigma[pro] == 0.4"))
grid$obs_sd_label = factor(grid$obs_sd, 
  labels = c("sigma[obs] == 0.2","sigma[obs] == 0.4","sigma[obs] == 0.8"))
grid$pro_CV_label = factor(grid$pro_CV, 
  labels = c("CV[pro] == 0.1","CV[pro] == 0.5","CV[pro] == 1"))
grid$obs_CV_label = factor(grid$obs_CV, 
  labels = c("CV[obs] == 0.1","CV[obs] == 0.5","CV[obs] == 1"))

# Figure 01 - observation error
pdf("plots/Figure_01_estimated_obs_error.pdf")
g1 = ggplot(dplyr::filter(grid, b_CV==1),
  aes(x=as.factor(obs_CV), y=sd_obs_est, group=as.factor(obs_CV))) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.65,outlier.shape = NA) +
  xlab(expression(prior~CV~sigma[obs])) + ylab(expression(Estimated~sigma[obs])) +
  geom_hline(aes(yintercept = obs_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label ~ pro_sd_label, scale="free_y",labeller = "label_parsed") + 
  theme_bw() + theme(strip.background = element_rect(color="black",fill="white"))
g1
dev.off()
# Figure 01 - same plot as Jpeg
jpeg("plots/Figure_01_estimated_obs_error.jpeg")
g1
dev.off()

pdf("plots/Figure_02_estimated_process_error.pdf")
g2 = ggplot(dplyr::filter(grid, b_CV==1),
  aes(x=as.factor(pro_CV), y=sd_pro_est, group=as.factor(pro_CV))) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.65,outlier.shape = NA) +
  xlab(expression(prior~CV~sigma[pro])) + ylab(expression(Estimated~sigma[pro])) +
  geom_hline(aes(yintercept = pro_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label~ pro_sd_label, scale="free_y",labeller = "label_parsed") + 
  theme_bw() + theme(strip.background = element_rect(color="black",fill="white"))
g2
dev.off()

jpeg("plots/Figure_02_estimated_process_error.jpeg")
g2
dev.off()

pdf("plots/Figure_03_obs_v_process_error.pdf")
g3 = ggplot(dplyr::filter(dplyr::filter(grid, b_CV==1), obs_sd == 0.2, pro_sd==0.2),
  aes(x=sd_obs_est, y=sd_pro_est)) +
  geom_point(col="darkblue",fill="darkblue",alpha=0.65,size=2) +
  facet_grid(obs_CV_label ~ pro_CV_label, scale="free",labeller = "label_parsed") +
  xlab(expression(Estimated~sigma[obs])) + ylab(expression(Estimated~sigma[pro])) +
  geom_vline(aes(xintercept = 0.2),col="red",alpha=0.3) +
  geom_hline(aes(yintercept = 0.2),col="red",alpha=0.3) + 
  theme_bw() + theme(strip.background = element_rect(color="black",fill="white"))
g3
dev.off()

jpeg("plots/Figure_03_obs_v_process_error.jpeg")
g3
dev.off()
