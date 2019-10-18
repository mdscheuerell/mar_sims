
library(ggplot2)
grid = readRDS("grid.rds")
post = readRDS("results/posterior_summaries.rds")

sd_obs_locs = grep("SD_obs",rownames(post))
grid$sd_obs_est = post$mean[sd_obs_locs]
sd_pro_locs = grep("SD_pro",rownames(post))
grid$sd_pro_est = post$mean[sd_pro_locs]

# make basic plots of 
pdf("plots/Estimated obs error.pdf")
ggplot(grid, 
  aes(x=as.factor(obs_CV), y=sd_obs_est, group=as.factor(obs_CV))) + 
  geom_boxplot() + xlab("prior CV of obs SD") + ylab("Estimated obs error") + 
  geom_hline(aes(yintercept = obs_sd),col="red") + 
  facet_grid(obs_sd~ pro_sd, scale="free")
dev.off()

pdf("plots/Estimated process error.pdf")
ggplot(grid, 
  aes(x=as.factor(pro_CV), y=sd_pro_est, group=as.factor(pro_CV))) + 
  geom_boxplot() + xlab("prior CV of process SD") + ylab("Estimated pro error") + 
  geom_hline(aes(yintercept = pro_sd),col="red") + 
  facet_grid(pro_sd~ obs_sd, scale="free")
dev.off()

pdf("plots/Obs v process error.pdf")
ggplot(dplyr::filter(grid, obs_sd == 0.2, pro_sd==0.2), 
  aes(x=sd_obs_est, y=sd_pro_est)) + 
  geom_point() + 
  facet_grid(obs_CV ~ pro_CV, scale="free") + 
  xlab("prior CV of obs SD") + ylab("prior CV of process SD")
dev.off()

grid$B11 = post$mean[which(substr(rownames(post),1,9)=="Bmat[1,1]")]
grid$B12 = post$mean[which(substr(rownames(post),1,9)=="Bmat[1,2]")]
grid$B21 = post$mean[which(substr(rownames(post),1,9)=="Bmat[2,1]")]
grid$B22 = post$mean[which(substr(rownames(post),1,9)=="Bmat[2,2]")]
grid$B23 = post$mean[which(substr(rownames(post),1,9)=="Bmat[2,3]")]
grid$B32 = post$mean[which(substr(rownames(post),1,9)=="Bmat[3,2]")]
grid$B33 = post$mean[which(substr(rownames(post),1,9)=="Bmat[3,3]")]
grid$B34 = post$mean[which(substr(rownames(post),1,9)=="Bmat[3,4]")]
grid$B43 = post$mean[which(substr(rownames(post),1,9)=="Bmat[4,3]")]
grid$B44 = post$mean[which(substr(rownames(post),1,9)=="Bmat[4,4]")]

# look at estimates of B[1,1] in same format as above
pdf("plots/Estimated B11.pdf")
ggplot(grid, 
  aes(x=as.factor(obs_CV), y=B11, group=as.factor(obs_CV))) + 
  geom_boxplot() + xlab("prior CV of obs SD") + ylab("Estimated B[1,1]") + 
  geom_hline(aes(yintercept = 0.5),col="red") + 
  facet_grid(obs_sd~ pro_sd, scale="free")
dev.off()

# similarly there's tradeoffs between B[1,1] and B[1,2]
pdf("plots/B11 v B12.pdf")
ggplot(dplyr::filter(grid, obs_sd == 0.2, pro_sd==0.2), 
  aes(x=B11, y=B12)) + 
  geom_point() + xlab("Estimated B[1,1]") + ylab("Estimated B[1,2]") + 
  facet_grid(obs_CV~ pro_CV, scale="free")
dev.off()

