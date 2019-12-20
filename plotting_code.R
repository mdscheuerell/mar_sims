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
grid$pro_sd_label = paste0("pro_sd=",grid$pro_sd)
grid$obs_sd_label = paste0("obs_sd=",grid$obs_sd)
grid$pro_CV_label = paste0("pro_CV=",grid$pro_CV)
grid$obs_CV_label = paste0("obs_CV=",grid$obs_CV)

# make basic plots of
pdf("plots/Estimated obs error.pdf")
ggplot(dplyr::filter(grid, b_CV==1),
  aes(x=as.factor(obs_CV), y=sd_obs_est, group=as.factor(obs_CV))) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  xlab("prior CV of obs SD") + ylab("Estimated obs error SD") +
  geom_hline(aes(yintercept = obs_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label~ pro_sd_label, scale="free")
dev.off()

pdf("plots/Estimated process error.pdf")
ggplot(dplyr::filter(grid, b_CV==1),
  aes(x=as.factor(pro_CV), y=sd_pro_est, group=as.factor(pro_CV))) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  xlab("prior CV of process SD") + ylab("Estimated pro error SD") +
  geom_hline(aes(yintercept = pro_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label~ pro_sd_label, scale="free")
dev.off()

pdf("plots/Obs v process error.pdf")
ggplot(dplyr::filter(dplyr::filter(grid, b_CV==1), obs_sd == 0.2, pro_sd==0.2),
  aes(x=sd_obs_est, y=sd_pro_est)) +
  geom_point() +
  facet_grid(obs_CV_label ~ pro_CV_label, scale="free") +
  xlab("Estimated obs error SD") + ylab("Estimated pro error SD") +
  geom_vline(aes(xintercept = 0.2),col="red",alpha=0.3) +
  geom_hline(aes(yintercept = 0.2),col="red",alpha=0.3)
dev.off()

# make similar plots for interactions
grid$iter = seq(1,nrow(grid))
post = dplyr::left_join(post, grid)

post$pro_sd_label = paste0("pro_sd=",post$pro_sd)
post$obs_sd_label = paste0("obs_sd=",post$obs_sd)
post$pro_CV_label = paste0("pro_CV=",post$pro_CV)
post$obs_CV_label = paste0("obs_CV=",post$obs_CV)

# Make plots of b as function of obs error CV
maxiter = max(post$iter)
b11 = ggplot(dplyr::filter(post, obs_sd==0.2, pro_sd==0.2,b_CV==1, pro_CV==0.1,
                     par %in% paste0("Bmat[1,1]",seq(1,maxiter))),
       aes(x=as.factor(obs_CV), y=mean,group=as.factor(obs_CV))) +
  geom_boxplot() + geom_hline(aes(yintercept = 0.5),col="red",alpha=0.3) +
  ylab("B[1,1] estimate") + xlab("Obs CV")
b12 = ggplot(dplyr::filter(post, obs_sd==0.2, pro_sd==0.2,b_CV==1, pro_CV==0.1,
                           par %in% paste0("Bmat[1,2]",seq(1,maxiter))),
             aes(x=as.factor(obs_CV), y=mean,group=as.factor(obs_CV))) +
  geom_boxplot() + geom_hline(aes(yintercept = 0.3),col="red",alpha=0.3) +
  ylab("estimate") + xlab("Obs CV")
b22 = ggplot(dplyr::filter(post, obs_sd==0.2, pro_sd==0.2,b_CV==1, pro_CV==0.1,
                           par %in% paste0("Bmat[2,2]",seq(1,maxiter))),
             aes(x=as.factor(obs_CV), y=mean,group=as.factor(obs_CV))) +
  geom_boxplot() + geom_hline(aes(yintercept = 0.6),col="red",alpha=0.3) +
  ylab("estimate") + xlab("Obs CV")
b21 = ggplot(dplyr::filter(post, obs_sd==0.2, pro_sd==0.2,b_CV==1, pro_CV==0.1,
                           par %in% paste0("Bmat[2,1]",seq(1,maxiter))),
             aes(x=as.factor(obs_CV), y=mean,group=as.factor(obs_CV))) +
  geom_boxplot() + geom_hline(aes(yintercept = -0.1),col="red",alpha=0.3) +
  ylab("estimate") + xlab("Obs CV")

pdf("Estimated_b_elements_as_function_of_obs_cv.pdf")
gridExtra::grid.arrange(b11,b12,b21,b22,nrow=2,ncol=2)
dev.off()

B0_lfc <- matrix(c(0.5, -0.1,  0.0,  0.0,
                   0.3,  0.6, -0.2,  0.0,
                   0.0,  0.2,  0.7, -0.3,
                   0.0,  0.0,  0.1,  0.8),4,4, byrow=TRUE)
post$shortpar = NA
post$true = NA
for(i in 1:4) {
  for(j in 1:4) {
    g = grep(paste0("Bmat[",i,",",j,"]"),post$par, fixed=TRUE)
    post$true[g] = B0_lfc[i,j]
    post$shortpar[g] = paste0("B[",i,",",j,"]")
  }
}

pdf("plots/Estimated_b_elements_bCV.pdf")
ggplot(dplyr::filter(post, obs_CV==1, pro_CV==1, !is.na(shortpar)),
       aes(x=as.factor(b_CV), y=mean,group=as.factor(b_CV))) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  ylab("Estimated B parameter") + xlab("prior CV on B parameter") +
  facet_wrap(~shortpar,scale="free_y") +
  geom_hline(aes(yintercept = true),col="red",alpha=0.3)
dev.off()

# make groups obs sd
pdf("plots/Estimated_b_elements_bCV_grouped.pdf")
post2 = dplyr::filter(post, obs_CV==1, pro_CV==1, !is.na(shortpar))
post2$group = NA
post2$group[which(post2$obs_sd==0.2 & post2$pro_sd == 0.1)] = "obs_sd = 0.2, pro_sd = 0.1"
post2$group[which(post2$obs_sd==0.2 & post2$pro_sd == 0.4)] = "obs_sd = 0.2, pro_sd = 0.4"
post2$group[which(post2$obs_sd==0.8 & post2$pro_sd == 0.1)] = "obs_sd = 0.8, pro_sd = 0.1"
post2$group[which(post2$obs_sd==0.8 & post2$pro_sd == 0.4)] = "obs_sd = 0.8, pro_sd = 0.4"
ggplot(dplyr::filter(post2,!is.na(group)),
       aes(x=as.factor(b_CV), y=mean,fill=group)) +
  geom_boxplot(alpha=0.4,outlier.shape = NA) +
  ylab("Estimated B parameter") + xlab("prior CV on B parameter") +
  facet_wrap(~shortpar,scale="free_y") +
  geom_hline(aes(yintercept = true),col="red",alpha=0.3)
dev.off()