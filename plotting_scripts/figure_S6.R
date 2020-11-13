library(dplyr)
library(ggplot2)
library(viridis)

## load grid of sim options
grid <- readRDS("grid.rds")
#new_grid = readRDS("grid_extended.rds")
#grid = rbind(grid, new_grid)
indx = which(grid$obs_sd==0.2 & grid$pro_sd==0.4)
indx2 = which(grid$obs_sd==0.4 & grid$pro_sd==0.2)
grid = grid[c(indx,indx2),]

grid = dplyr::filter(grid, pro_CV==1, obs_CV==1, b_CV==1)

post = readRDS("results/posterior_summaries_2survey.rds")
post$par = rownames(post)

#grid = dplyr::filter(grid, iter<=max(post$iter))

sd_obs_locs = grep("SD_obs",rownames(post))
grid$sd_obs_est = post$mean[sd_obs_locs]
sd_pro_locs = grep("SD_pro",rownames(post))
grid$sd_pro_est = post$mean[sd_pro_locs]


# labels
grid$pro_sd_label = factor(grid$pro_sd, 
  labels = c("sigma[pro] == 0.2","sigma[pro] == 0.4"))
grid$obs_sd_label = factor(grid$obs_sd, 
  labels = c("sigma[obs] == 0.2","sigma[obs] == 0.4"))
grid$pro_CV_label = factor(grid$pro_CV, 
  labels = c("CV[pro] == 0.1","CV[pro] == 0.5","CV[pro] == 1"))
grid$obs_CV_label = factor(grid$obs_CV, 
  labels = c("CV[obs] == 0.1","CV[obs] == 0.5","CV[obs] == 1"))


# make similar plots for interactions
#grid$iter = seq(1,nrow(grid))
post = dplyr::left_join(post, grid)

# Make plots of b as function of obs error CV
maxiter = max(post$iter)

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


# Just show the two contrasting scenariors
post2 = dplyr::filter(post, 
  obs_CV==1, pro_CV==1,b_CV==1,
  !is.na(shortpar))

post2$group = as.factor(paste(post2$obs_sd, post2$pro_sd))
post2$group_label = factor(post2$group, 
  labels = c("sigma[obs] == 0.2, sigma[pro] == 0.4",
    "sigma[obs] == 0.4, sigma[pro] == 0.2"))


pdf("plots/Figure_S6.pdf")

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
dev.off()

jpeg("plots/Figure_S6.jpeg")
g1
dev.off()

