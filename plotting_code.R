library(dplyr)
library(ggplot2)
grid = readRDS("grid.rds")
post = rbind(readRDS("results/posterior_summaries_t_1.rds"),
             readRDS("results/posterior_summaries_t_2.rds"),
             readRDS("results/posterior_summaries_t_3.rds"))
post$par = rownames(post)

# these added later
newgrid = readRDS("grid_extended.rds")
newpost = readRDS("results/posterior_summaries_1_extended.rds")
newpost$par = rownames(newpost)

sd_obs_locs = grep("SD_obs",rownames(post))
grid$sd_obs_est = post$mean[sd_obs_locs]
sd_pro_locs = grep("SD_pro",rownames(post))
grid$sd_pro_est = post$mean[sd_pro_locs]

sd_obs_locs = grep("SD_obs",rownames(newpost))
newgrid$sd_obs_est = newpost$mean[sd_obs_locs]
sd_pro_locs = grep("SD_pro",rownames(newpost))
newgrid$sd_pro_est = newpost$mean[sd_pro_locs]

# labels
grid$pro_sd_label = factor(grid$pro_sd, 
  labels = c("sigma[pro] == 0.1","sigma[pro] == 0.2","sigma[pro] == 0.4"))
grid$obs_sd_label = factor(grid$obs_sd, 
  labels = c("sigma[obs] == 0.2","sigma[obs] == 0.4","sigma[obs] == 0.8"))
grid$pro_CV_label = factor(grid$pro_CV, 
  labels = c("CV[pro] == 0.1","CV[pro] == 0.5","CV[pro] == 1"))
grid$obs_CV_label = factor(grid$obs_CV, 
  labels = c("CV[obs] == 0.1","CV[obs] == 0.5","CV[obs] == 1"))

# labels
newgrid$pro_sd_label = factor(newgrid$pro_sd, 
  labels = c("sigma[pro] == 0.2"))
newgrid$obs_sd_label = factor(newgrid$obs_sd, 
  labels = c("sigma[obs] == 0.4"))
newgrid$pro_CV_label = factor(newgrid$pro_CV, 
  labels = c("CV[pro] == 0.1","CV[pro] == 0.5","CV[pro] == 1"))
newgrid$obs_CV_label = factor(newgrid$obs_CV, 
  labels = c("CV[obs] == 0.1","CV[obs] == 0.5","CV[obs] == 1"))

# make basic plots of
pdf("plots/Estimated obs error.pdf")
ggplot(dplyr::filter(grid, b_CV==1),
  aes(x=as.factor(obs_CV), y=sd_obs_est, group=as.factor(obs_CV))) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.65,outlier.shape = NA) +
  xlab(expression(prior~CV~sigma[obs])) + ylab(expression(Estimated~sigma[obs])) +
  geom_hline(aes(yintercept = obs_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label ~ pro_sd_label, scale="free_y",labeller = "label_parsed") + 
  theme_bw() + theme(strip.background = element_rect(color="black",fill="white"))
dev.off()

pdf("plots/Estimated process error.pdf")
ggplot(dplyr::filter(grid, b_CV==1),
  aes(x=as.factor(pro_CV), y=sd_pro_est, group=as.factor(pro_CV))) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.65,outlier.shape = NA) +
  xlab(expression(prior~CV~sigma[pro])) + ylab(expression(Estimated~sigma[pro])) +
  geom_hline(aes(yintercept = pro_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label~ pro_sd_label, scale="free_y",labeller = "label_parsed") + 
  theme_bw() + theme(strip.background = element_rect(color="black",fill="white"))
dev.off()

pdf("plots/Obs v process error.pdf")
ggplot(dplyr::filter(dplyr::filter(grid, b_CV==1), obs_sd == 0.2, pro_sd==0.2),
  aes(x=sd_obs_est, y=sd_pro_est)) +
  geom_point(col="darkblue",fill="darkblue",alpha=0.65,size=2) +
  facet_grid(obs_CV_label ~ pro_CV_label, scale="free",labeller = "label_parsed") +
  xlab(expression(Estimated~sigma[obs])) + ylab(expression(Estimated~sigma[pro])) +
  geom_vline(aes(xintercept = 0.2),col="red",alpha=0.3) +
  geom_hline(aes(yintercept = 0.2),col="red",alpha=0.3) + 
  theme_bw() + theme(strip.background = element_rect(color="black",fill="white")) + 
  
dev.off()

# make similar plots for interactions
grid$iter = seq(1,nrow(grid))
post = dplyr::left_join(post, grid)
newgrid$iter = seq(1,nrow(newgrid))
newpost = dplyr::left_join(newpost, newgrid)
#post$pro_sd_label = paste0("pro_sd=",post$pro_sd)
#post$obs_sd_label = paste0("obs_sd=",post$obs_sd)
#post$pro_CV_label = paste0("pro_CV=",post$pro_CV)
#post$obs_CV_label = paste0("obs_CV=",post$obs_CV)

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
newpost$shortpar = NA
newpost$true = NA
for(i in 1:4) {
  for(j in 1:4) {
    g = grep(paste0("Bmat[",i,",",j,"]"),post$par, fixed=TRUE)
    post$true[g] = B0_lfc[i,j]
    post$shortpar[g] = paste0("B[",i,",",j,"]")
    g = grep(paste0("Bmat[",i,",",j,"]"),newpost$par, fixed=TRUE)
    newpost$true[g] = B0_lfc[i,j]
    newpost$shortpar[g] = paste0("B[",i,",",j,"]")    
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

# Same plot as above, but ungrouped
pdf("plots/Estimated_b_elements_bCV_ungrouped.pdf")
post2 = dplyr::filter(post, obs_CV==1, pro_CV==1, !is.na(shortpar))
post2$group = NA
post2$group[which(post2$obs_sd==0.2 & post2$pro_sd == 0.1)] = "obs_sd = 0.2, pro_sd = 0.1"
post2$group[which(post2$obs_sd==0.2 & post2$pro_sd == 0.4)] = "obs_sd = 0.2, pro_sd = 0.4"
post2$group[which(post2$obs_sd==0.8 & post2$pro_sd == 0.1)] = "obs_sd = 0.8, pro_sd = 0.1"
post2$group[which(post2$obs_sd==0.8 & post2$pro_sd == 0.4)] = "obs_sd = 0.8, pro_sd = 0.4"

ggplot(dplyr::filter(post2,!is.na(group),post2$obs_sd==0.2,post2$pro_sd == 0.1),
  aes(x=as.factor(b_CV), y=mean)) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  ylab("Estimated B parameter") + xlab("prior CV on B parameter") +
  facet_wrap(~shortpar,scale="free_y") +
  geom_hline(aes(yintercept = true),col="red",alpha=0.3) + 
  ggtitle("obs_sd = 0.2, pro_sd = 0.1")

ggplot(dplyr::filter(post2,!is.na(group),post2$obs_sd==0.2,post2$pro_sd == 0.4),
  aes(x=as.factor(b_CV), y=mean)) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  ylab("Estimated B parameter") + xlab("prior CV on B parameter") +
  facet_wrap(~shortpar,scale="free_y") +
  geom_hline(aes(yintercept = true),col="red",alpha=0.3) + 
  ggtitle("obs_sd = 0.2, pro_sd = 0.4")

ggplot(dplyr::filter(post2,!is.na(group),post2$obs_sd==0.8,post2$pro_sd == 0.1),
  aes(x=as.factor(b_CV), y=mean)) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  ylab("Estimated B parameter") + xlab("prior CV on B parameter") +
  facet_wrap(~shortpar,scale="free_y") +
  geom_hline(aes(yintercept = true),col="red",alpha=0.3) + 
  ggtitle("obs_sd = 0.8, pro_sd = 0.1")

ggplot(dplyr::filter(post2,!is.na(group),post2$obs_sd==0.8,post2$pro_sd == 0.4),
  aes(x=as.factor(b_CV), y=mean)) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  ylab("Estimated B parameter") + xlab("prior CV on B parameter") +
  facet_wrap(~shortpar,scale="free_y") +
  geom_hline(aes(yintercept = true),col="red",alpha=0.3) + 
  ggtitle("obs_sd = 0.8, pro_sd = 0.4")
dev.off()


# Just show the two contrasting scenariors
post2 = dplyr::filter(post, 
  obs_CV==1, pro_CV==1, 
  !is.na(shortpar), 
  obs_sd == 0.2, pro_sd == 0.4)

# bring in new data
newpost2 = dplyr::filter(newpost, 
  obs_CV==1, pro_CV==1,
  !is.na(shortpar))

post2 = rbind(post2, newpost2)
post2$group = as.factor(paste(post2$obs_sd, post2$pro_sd))
post2$group_label = factor(post2$group, 
  labels = c("sigma[obs] == 0.2, sigma[pro] == 0.4",
    "sigma[obs] == 0.4, sigma[pro] == 0.2"))


pdf("plots/Estimated_b_elements_2scenarios.pdf")

post2$mean[which(post2$mean==0)] = NA
post2$true[which(is.na(post2$mean))] = NA
library(viridis)

ggplot(post2,
  aes(x = group, y=mean, group=group, col=group,fill=group)) +
  geom_boxplot(alpha=0.4,outlier.shape = NA) + 
  ylab("Estimated B parameter") + 
  xlab("Scenario") +
  facet_wrap(~shortpar,scale="free_y") +
  geom_hline(aes(yintercept = true),col="red",alpha=0.3) + 
  scale_fill_viridis(discrete = TRUE, end=0.7) + 
  scale_color_viridis(discrete = TRUE, end=0.7) + 
  theme_bw() + 
  theme(legend.position='none',strip.background = element_rect(color="black",fill="white")) + 
scale_x_discrete(labels = c('0.2 0.4' = expression(paste(sigma[obs],"= 0.2")),
  '0.4 0.2'   = expression(paste(sigma[obs],"= 0.4"))))
dev.off()


# Esrtimates of MLEs vs Posteriors

post = rbind(readRDS("results/posterior_summaries_t_1.rds"),
  readRDS("results/posterior_summaries_t_2.rds"),
  readRDS("results/posterior_summaries_t_3.rds"))
post$par = rownames(post)

grid = readRDS("marss_pars.rds")
grid$pro_sd_label = paste0("pro_sd=",grid$pro_sd)
grid$obs_sd_label = paste0("obs_sd=",grid$obs_sd)
grid$pro_CV_label = paste0("pro_CV=",grid$pro_CV)
grid$obs_CV_label = paste0("obs_CV=",grid$obs_CV)

grid$sd_obs_est = sqrt(grid$R)
grid$sd_pro_est = sqrt(grid$Q)
  

# make basic plots of
pdf("plots/Estimated obs error_marss.pdf")
ggplot(dplyr::filter(grid, b_CV==1),
       aes(y=sd_obs_est)) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  ylab("Estimated obs error SD") +
  geom_hline(aes(yintercept = obs_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label~ pro_sd_label, scale="free")
dev.off()

pdf("plots/Estimated pro error_marss.pdf")
ggplot(dplyr::filter(grid, b_CV==1),
  aes(y=sd_pro_est)) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  ylab("Estimated pro error SD") +
  geom_hline(aes(yintercept = pro_sd),col="red",alpha=0.3) +
  facet_grid(obs_sd_label~ pro_sd_label, scale="free")
dev.off()

pdf("plots/Obs v process error_marss.pdf")
ggplot(dplyr::filter(dplyr::filter(grid, b_CV==1)),
  aes(x=sd_obs_est, y=sd_pro_est)) +
  geom_point() +
  facet_grid(obs_sd_label ~ pro_sd_label, scale="free") +
  xlab("Estimated obs error SD") + ylab("Estimated pro error SD") +
  geom_vline(aes(xintercept = obs_sd),col="red",alpha=0.3) +
  geom_hline(aes(yintercept = pro_sd),col="red",alpha=0.3)
dev.off()


# make similar plots for interactions
grid$iter = seq(1,nrow(grid))

# start with upper right corner or obs_erorr plot, pro_sd = 0.4, obs_sd = 0.2
post2 = post[grep("SD_obs",post$par),]
post2 = dplyr::left_join(grid,post2)

post2$obs_CV = as.factor(post2$obs_CV)
post2$pro_CV_label = paste("pro_sd = ",post2$pro_sd)
pdf("marss_v_bayes_sigma_obs.pdf")
ggplot(filter(post2,b_CV==1, pro_sd==0.4, obs_sd==0.2), aes(mean, sd_obs_est, col=obs_CV)) + 
  geom_point() + 
  facet_wrap(~pro_CV_label) + 
  geom_abline(intercept=0,slope=1) + 
  xlab("Posterior mean, sigma_obs") + 
  ylab("MLE, sigma_obs") #+ 
  #ggtitle("MLE tends to lead to higher estimates of sigma_obs when pro_sd > obs_sd")
dev.off()

# start with upper right corner or obs_erorr plot, pro_sd = 0.4, obs_sd = 0.2
post2 = post[grep("SD_pro",post$par),]
post2 = dplyr::left_join(grid,post2)

post2$pro_CV = as.factor(post2$pro_CV)
post2$obs_CV_label = paste("obs_sd = ",post2$obs_sd)
pdf("marss_v_bayes_sigma_pro.pdf")
ggplot(filter(post2,b_CV==1, pro_sd==0.4, obs_sd==0.2), aes(mean, sd_pro_est, col=pro_CV)) + 
  geom_point() + 
  facet_wrap(~obs_CV_label) + 
  geom_abline(intercept=0,slope=1) + 
  xlab("Posterior mean, sigma_obs") + 
  ylab("MLE, sigma_pro") + 
  ggtitle("MLE tends to lead to lower estimates of sigma_obs when pro_sd > obs_sd")
dev.off()


pdf("marss_v_bayes_sigma_obs.pdf")

post2 = post[which(startsWith(post$par, "Bmat[1,1]")==TRUE),]
post2 = dplyr::left_join(grid,post2)
post2$b_CV = as.factor(post2$b_CV)
post2$b_CV_label = paste("pro_sd = ",post2$pro_sd)

post2=filter(post2,pro_sd==0.4, obs_sd==0.2,pro_CV==1,obs_CV==1)
b11 = ggplot(post2, aes(mean, b11, col=b_CV)) + 
  geom_point() + 
  geom_abline(intercept=0,slope=1) + 
  xlab("Posterior mean, B[1,1]") + 
  ylab("MLE, B[1,1]")

post2 = post[which(startsWith(post$par, "Bmat[1,2]")==TRUE),]
post2 = dplyr::left_join(grid,post2)
post2$b_CV = as.factor(post2$b_CV)
post2$b_CV_label = paste("pro_sd = ",post2$pro_sd)

post2=filter(post2,pro_sd==0.4, obs_sd==0.2,pro_CV==1,obs_CV==1)
b12 = ggplot(post2, aes(mean, b12, col=b_CV)) + 
  geom_point() + 
  geom_abline(intercept=0,slope=1) + 
  xlab("Posterior mean, B[1,2]") + 
  ylab("MLE, B[1,2]")

post2 = post[which(startsWith(post$par, "Bmat[2,1]")==TRUE),]
post2 = dplyr::left_join(grid,post2)
post2$b_CV = as.factor(post2$b_CV)
post2$b_CV_label = paste("pro_sd = ",post2$pro_sd)

post2=filter(post2,pro_sd==0.4, obs_sd==0.2,pro_CV==1,obs_CV==1)
b21 = ggplot(post2, aes(mean, b21, col=b_CV)) + 
  geom_point() + 
  geom_abline(intercept=0,slope=1) + 
  xlab("Posterior mean, B[2,1]") + 
  ylab("MLE, B[2,1]")

post2 = post[which(startsWith(post$par, "Bmat[2,2]")==TRUE),]
post2 = dplyr::left_join(grid,post2)
post2$b_CV = as.factor(post2$b_CV)
post2$b_CV_label = paste("pro_sd = ",post2$pro_sd)

post2=filter(post2,pro_sd==0.4, obs_sd==0.2,pro_CV==1,obs_CV==1)
b22 = ggplot(post2, aes(mean, b22, col=b_CV)) + 
  geom_point() + 
  geom_abline(intercept=0,slope=1) + 
  xlab("Posterior mean, B[2,2]") + 
  ylab("MLE, B[2,2]")

gridExtra::grid.arrange(b11,b12,b21,b22,nrow=2,ncol=2)
dev.off()