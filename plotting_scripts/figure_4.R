library(dplyr)
library(ggplot2)
library(viridis)

grid = readRDS("grid.rds")
post = rbind(readRDS("results/posterior_summaries_t_1.rds"),
             readRDS("results/posterior_summaries_t_2.rds"),
             readRDS("results/posterior_summaries_t_3.rds"))
post$par = rownames(post)

# these added for the contrasting scenarios
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

# make similar plots for interactions
grid$iter = seq(1,nrow(grid))
post = dplyr::left_join(post, grid)
newgrid$iter = seq(1,nrow(newgrid))
newpost = dplyr::left_join(newpost, newgrid)

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

# Just show the two contrasting scenariors
post2 = dplyr::filter(post, 
  obs_CV==1, pro_CV==1, b_CV==1,
  !is.na(shortpar), 
  obs_sd == 0.2, pro_sd == 0.4)

# bring in new data
newpost2 = dplyr::filter(newpost, 
  obs_CV==1, pro_CV==1,b_CV==1,
  !is.na(shortpar))

post2 = rbind(post2, newpost2)
post2$group = as.factor(paste(post2$obs_sd, post2$pro_sd))
post2$group_label = factor(post2$group, 
  labels = c("sigma[obs] == 0.2, sigma[pro] == 0.4",
    "sigma[obs] == 0.4, sigma[pro] == 0.2"))


pdf("plots/Figure_4.pdf")

post2$mean[which(post2$mean==0)] = NA
post2$true[which(is.na(post2$mean))] = NA


g1 = ggplot(post2,
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
  scale_x_discrete(labels = c('0.2 0.4' = expression(atop(paste(sigma[obs],"= 0.2"), paste(sigma[pro],"= 0.4"))),
    '0.4 0.2'   = expression(atop(paste(sigma[obs],"= 0.4"),paste(sigma[pro],"= 0.2")))))
g1
dev.off()

jpeg("plots/Figure_4.jpeg")
g1
dev.off()
# Esrtimates of MLEs vs Posteriors
