library(dplyr)
library(ggplot2)
library(viridis)

# filter the grid to the 2 scenarios of interest
grid = readRDS("grid.rds")
indx = which(grid$obs_sd==0.2 & grid$pro_sd==0.4)
indx2 = which(grid$obs_sd==0.4 & grid$pro_sd==0.2)
grid = grid[c(indx,indx2),]
grid = dplyr::filter(grid, pro_CV==1, obs_CV==1, b_CV==1)

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
x1 = dplyr::filter(post, seed %in% keep$seed,
                     seed %in% grid$seed)
x1$Surveys = 1
x1 = x1[grep("Bmat",x1$par),]
x1 = dplyr::left_join(x1, grid)

# Add in the 2 scenario data
x = readRDS("results/posterior_summaries_2surveys.rds")
x = x[grep("Bmat",x$par),]
x2 = dplyr::left_join(x, grid)
x2$Surveys = 2
keep = dplyr::group_by(x2, seed) %>%
  dplyr::summarize(max_rhat = max(Rhat,na.rm=T)) %>%
  dplyr::filter(max_rhat <= 1.1)
x2 = x2[which(x2$seed %in% keep$seed),]

# Add in the 4 scenario data
x = readRDS("results/posterior_summaries_4surveys.rds")
x = x[grep("Bmat",x$par),]
x4 = dplyr::left_join(x, grid)
x4$Surveys = 4
keep = dplyr::group_by(x4, seed) %>%
  dplyr::summarize(max_rhat = max(Rhat,na.rm=T)) %>%
  dplyr::filter(max_rhat <= 1.1)
x4 = x4[which(x4$seed %in% keep$seed),]

# Join the 3 together
x = rbind(x1, x2, x4)

# Add names for scenarios
x$Scenario = paste0(x$obs_sd, ", ", x$pro_sd)

# Add true values for b matrix
x$true = NA
x$true[which(x$par=="Bmat[1,1]")] = 0.5
x$true[which(x$par=="Bmat[1,2]")] = -0.1
x$true[which(x$par=="Bmat[2,1]")] = 0.3
x$true[which(x$par=="Bmat[2,2]")] = 0.6
x$true[which(x$par=="Bmat[2,3]")] = -0.2
x$true[which(x$par=="Bmat[3,2]")] = 0.2
x$true[which(x$par=="Bmat[3,3]")] = 0.7
x$true[which(x$par=="Bmat[3,4]")] = -0.3
x$true[which(x$par=="Bmat[4,3]")] = 0.1
x$true[which(x$par=="Bmat[4,4]")] = 0.8

x$mean[which(x$par %in% c("Bmat[1,3]","Bmat[1,4]","Bmat[2,4]","Bmat[3,1]","Bmat[4,1]","Bmat[4,2]"))] = NA

x$median = x[,"50%"]
x$Surveys = as.factor(x$Surveys)

x$group = as.factor(paste(x$Surveys,x$Scenario))
x$new_group = as.numeric(x$group)
x$new_group[which(as.numeric(x$group)%in%c(1,3,5))] = x$new_group[which(as.numeric(x$group)%in%c(1,3,5))] + 0.25
x$new_group[which(as.numeric(x$group)%in%c(2,4,6))] = x$new_group[which(as.numeric(x$group)%in%c(2,4,6))] - 0.25
x$new_group = x$new_group - 0.25 + 0.005*x$iter

g1 = ggplot(x[x$iter < 95,], aes(new_group,mean, group=group, col=Scenario, fill=Scenario)) +
  #geom_boxplot(alpha=0.3,outlier.shape = NA) +
  #geom_point(position=position_dodge(width=0.75),aes(group=Scenario),alpha=0.1)+
  geom_linerange(aes(ymin=`25%`,ymax=`75%`),alpha=0.1,outlier.shape = NA) +
  geom_boxplot(alpha=0.3,fill=NA,outlier.shape = NA) +
  ylab("Estimated B parameter") +
  xlab("Surveys")+
  facet_wrap(~par,scale="free_y") +
  geom_hline(aes(yintercept = true),col="red",alpha=0.3) +
  scale_fill_viridis(discrete = TRUE, end=0.7,labels = expression(atop(paste(sigma[obs],"= 0.2"), paste(sigma[pro],"= 0.4")), atop(paste(sigma[obs],"= 0.4"), paste(sigma[pro],"= 0.2")) )) +
  scale_color_viridis(discrete = TRUE, end=0.7,labels = expression(atop(paste(sigma[obs],"= 0.2"), paste(sigma[pro],"= 0.4")), atop(paste(sigma[obs],"= 0.4"), paste(sigma[pro],"= 0.2")) )) +
  scale_x_continuous(breaks=c(1.5,3.5,5.5), labels = c("1","2","4")) +  
  theme_bw() +
  theme(strip.background = element_rect(color="black",fill="white"))


pdf("plots/Figure_3.pdf")
g1
dev.off()

jpeg("plots/Figure_3.jpeg")
g1
dev.off()


###############
# Make coverage plot for these 3 scenarios
x$coverage = 0
indx = which((x$true < x$`90%`) & (x$true > x$`10%`))
x$coverage[indx] = 1

g = dplyr::group_by(x, Surveys, Scenario, par) %>%
      dplyr::summarise(cover = mean(coverage))

g$cover[which(g$par %in% c("Bmat[1,3]","Bmat[1,4]","Bmat[2,4]","Bmat[3,1]","Bmat[4,1]","Bmat[4,2]"))] = NA

p <- ggplot(g, aes(Surveys,cover, col=Scenario)) +
  geom_point(size=2.5)+
  ylab("Coverage of posterior CIs") +
  xlab("Scenario") +
  facet_wrap(~par) +
  geom_hline(aes(yintercept=0.8),col="red",alpha=0.7) +
  scale_color_viridis(discrete = TRUE, end=0.7,labels = expression(atop(paste(sigma[obs],"= 0.2"), paste(sigma[pro],"= 0.4")), atop(paste(sigma[obs],"= 0.4"), paste(sigma[pro],"= 0.2")) )) +
  theme_bw() +
  theme(strip.background = element_rect(color="black",fill="white"))

pdf("plots/Figure_S7_coverage.pdf")
p
dev.off()

jpeg("plots/Figure_S7_coverage.jpeg")
p
dev.off()

############## Also show the standard deviation
g = dplyr::group_by(x, Surveys, Scenario, par) %>%
  dplyr::summarise(mean_sd = mean(sd))

g$mean_sd[which(g$par %in% c("Bmat[1,3]","Bmat[1,4]","Bmat[2,4]","Bmat[3,1]","Bmat[4,1]","Bmat[4,2]"))] = NA

p <- ggplot(g, aes(Surveys,mean_sd, col=Scenario)) +
  geom_point(size=2.5)+
  ylab("Standard deviation of parameter") +
  xlab("Scenario") +
  facet_wrap(~par, scale="free_y") +
  scale_color_viridis(discrete = TRUE, end=0.7,labels = expression(atop(paste(sigma[obs],"= 0.2"), paste(sigma[pro],"= 0.4")), atop(paste(sigma[obs],"= 0.4"), paste(sigma[pro],"= 0.2")) )) +
  theme_bw() +
  theme(strip.background = element_rect(color="black",fill="white"))

pdf("plots/Figure_S6_sd.pdf")
p
dev.off()

jpeg("plots/Figure_S6_sd.jpeg")
p
dev.off()
