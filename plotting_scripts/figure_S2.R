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

#mean(c(0.431,0.464,0.534,0.6)/seq(0.5,0.8,length.out=4))
dplyr::filter(post, obs_CV==1, pro_CV==1, !is.na(shortpar)) %>%
  dplyr::group_by(b_CV,shortpar) %>% 
  dplyr::summarize(m = mean(mean)) %>%
  dplyr::filter(b_CV==1, shortpar %in% c("B[1,1]","B[2,2]","B[3,3]","B[4,4]"))
  
pdf("plots/Figure_S2_estimated_b_elements_bCV.pdf")
g5 = ggplot(dplyr::filter(post, obs_CV==1, pro_CV==1, !is.na(shortpar)),
       aes(x=as.factor(b_CV), y=mean,group=as.factor(b_CV))) +
  geom_boxplot(col="darkblue",fill="darkblue",alpha=0.4,outlier.shape = NA) +
  ylab("Estimated B parameter") + xlab("prior CV on B parameter") +
  facet_wrap(~shortpar,scale="free_y") +
  geom_hline(aes(yintercept = true),col="red",alpha=0.3) + 
  theme_bw() + 
  theme(legend.position='none',strip.background = element_rect(color="black",fill="white"))
g5
dev.off()

jpeg("plots/Figure_S2_estimated_b_elements_bCV.jpeg")
g5
dev.off()
