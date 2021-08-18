library(dplyr)
library(ggplot2)
library(viridis)
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


# labels
grid$pro_sd_label = factor(grid$pro_sd,
  labels = c("sigma[pro] == 0.1","sigma[pro] == 0.2","sigma[pro] == 0.4"))
grid$obs_sd_label = factor(grid$obs_sd,
  labels = c("sigma[obs] == 0.2","sigma[obs] == 0.4","sigma[obs] == 0.8"))
grid$pro_CV_label = factor(grid$pro_CV,
  labels = c("CV[pro] == 0.1","CV[pro] == 0.5","CV[pro] == 1"))
grid$obs_CV_label = factor(grid$obs_CV,
  labels = c("CV[obs] == 0.1","CV[obs] == 0.5","CV[obs] == 1"))

# make similar plots for interactions
#grid$iter = seq(1,nrow(grid))
post = dplyr::left_join(post, grid)

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

post_summary = dplyr::filter(post,
  obs_CV==1,
  pro_CV==1,
  b_CV==1,
  !is.na(shortpar))

post_summary$group = as.factor(paste(post_summary$obs_sd, post_summary$pro_sd))
post_summary = dplyr::filter(post_summary,
  group %in% c("0.2 0.4","0.4 0.2"))
post_summary$group = as.factor(as.character(post_summary$group))

post_summary$group_label = factor(post_summary$group,
  labels = c("sigma[obs] == 0.2, sigma[pro] == 0.4",
    "sigma[obs] == 0.4, sigma[pro] == 0.2"))

pdf("plots/Figure_2b.pdf")

post_summary$mean[which(post_summary$mean==0)] = NA
post_summary$true[which(is.na(post_summary$mean))] = NA

post_summary$new_group = as.numeric(as.factor(post_summary$group))
post_summary$new_group = post_summary$new_group - 0.25 + 0.005*post_summary$iter

g1 = ggplot(post_summary[which(post_summary$iter < 95),],
            aes(x = new_group, y=mean, group=group, col=group,fill=group)) +
  geom_linerange(aes(ymin=`25%`,ymax=`75%`),alpha=0.1,outlier.shape = NA) +
  geom_boxplot(alpha=0.3,outlier.shape = NA, fill=NA)+
  #geom_sina(size=1, alpha=0.35) +
  #geom_violin(alpha=0.2,outlier.shape = NA,draw_quantiles = c(0.25, 0.5, 0.75)) +
  ylab("Estimated B parameter") +
  xlab("Scenario") +
  facet_wrap(~shortpar,scale="free_y") +
  geom_hline(aes(yintercept = true),col="red",alpha=0.3) +
  scale_fill_viridis(discrete = TRUE, end=0.7) +
  scale_color_viridis(discrete = TRUE, end=0.7) +
  theme_bw() +
  theme(legend.position='none',strip.background = element_rect(color="black",fill="white")) +
  scale_x_continuous(breaks=c(1,2), 
                     labels = c(expression(atop(paste(sigma[obs],"= 0.2"), paste(sigma[pro],"= 0.4"))),
                                expression(atop(paste(sigma[obs],"= 0.4"), paste(sigma[pro],"= 0.2"))))) 
g1
dev.off()

jpeg("plots/Figure_2b.jpeg")
g1
dev.off()
# Esrtimates of MLEs vs Posteriors
