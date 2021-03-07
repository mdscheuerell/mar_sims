library(dplyr)
library(ggplot2)
library(ggforce)

grid = readRDS("grid.rds")

# load in batched results
for(i in 1:9) {
  test = readRDS(paste0("results/posterior_summaries_",i,".rds"))
  if(i == 1) {
    post = test
  } else {
    post = rbind(post, test)
  }
}
# filter runs based on convergence --
keep = dplyr::group_by(post, seed) %>%
  dplyr::summarize(max_rhat = max(Rhat,na.rm=T)) %>%
  dplyr::filter(max_rhat <= 1.1)
post = dplyr::filter(post, seed %in% keep$seed)

sd_obs_locs = grep("SD_obs",post$par)
obs_sd_dat = post[sd_obs_locs,] %>%
  dplyr::select(mean,seed) %>%
  dplyr::rename(sd_obs_est = mean)
grid = dplyr::left_join(grid, obs_sd_dat)

sd_pro_locs = grep("SD_pro",post$par)
pro_sd_dat = post[sd_pro_locs,] %>%
  dplyr::select(mean,seed) %>%
  dplyr::rename(sd_pro_est = mean)
grid = dplyr::left_join(grid, pro_sd_dat)

# labels
grid$pro_sd_label = factor(grid$pro_sd,
  labels = c("sigma[pro] == 0.1","sigma[pro] == 0.2","sigma[pro] == 0.4"))
grid$obs_sd_label = factor(grid$obs_sd,
  labels = c("sigma[obs] == 0.2","sigma[obs] == 0.4","sigma[obs] == 0.8"))
grid$pro_CV_label = factor(grid$pro_CV,
  labels = c("CV[pro] == 0.1","CV[pro] == 0.5","CV[pro] == 1"))
grid$obs_CV_label = factor(grid$obs_CV,
  labels = c("CV[obs] == 0.1","CV[obs] == 0.5","CV[obs] == 1"))

# Figure 04 - observation error
pdf("plots/Figure_04_estimated_obs_error.pdf")
g1 = ggplot(dplyr::filter(grid, b_CV==1),
  aes(x=as.factor(obs_CV), y=sd_obs_est, group=as.factor(obs_CV))) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.3,outlier.shape = NA) +
  geom_point(col="darkblue", alpha=0.1) +
  #geom_violin(col="darkblue",fill="darkblue",alpha=0.2,outlier.shape = NA,draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_sina(col = "darkblue", size=1, alpha=0.35) +
  xlab(expression(prior~CV~sigma[obs])) + ylab(expression(Estimated~sigma[obs])) +
  geom_hline(aes(yintercept = obs_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label ~ pro_sd_label, scale="free_y",labeller = "label_parsed") +
  theme_bw() + theme(strip.background = element_rect(color="black",fill="white"))
g1
dev.off()
# Figure 04 - same plot as Jpeg
jpeg("plots/Figure_04_estimated_obs_error.jpeg")
g1
dev.off()

pdf("plots/Figure_05_estimated_process_error.pdf")
g2 = ggplot(dplyr::filter(grid, b_CV==1),
  aes(x=as.factor(pro_CV), y=sd_pro_est, group=as.factor(pro_CV))) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  geom_point(col="darkblue", alpha=0.1) +
  #geom_violin(col="darkblue",fill="darkblue",alpha=0.2,outlier.shape = NA,draw_quantiles = c(0.25, 0.5, 0.75)) +
  #geom_sina(col = "darkblue", size=1, alpha=0.35) +
  xlab(expression(prior~CV~sigma[pro])) + ylab(expression(Estimated~sigma[pro])) +
  geom_hline(aes(yintercept = pro_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label~ pro_sd_label, scale="free_y",labeller = "label_parsed") +
  theme_bw() + theme(strip.background = element_rect(color="black",fill="white"))
g2
dev.off()

jpeg("plots/Figure_05_estimated_process_error.jpeg")
g2
dev.off()

pdf("plots/Figure_06_obs_v_process_error.pdf")
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

jpeg("plots/Figure_06_obs_v_process_error.jpeg")
g3
dev.off()
