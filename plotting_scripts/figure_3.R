library(dplyr)
library(ggplot2)
library(viridis)

# filter the grid to the 2 scenarios of interest
grid = readRDS("grid.rds")
indx = which(grid$obs_sd==0.2 & grid$pro_sd==0.4)
indx2 = which(grid$obs_sd==0.4 & grid$pro_sd==0.2)
grid = grid[c(indx,indx2),]
grid = dplyr::filter(grid, pro_CV==1, obs_CV==1, b_CV==1)
grid$iter = seq(1,nrow(grid))

x = readRDS("results/posterior_summaries_2survey.rds")
x$par = rownames(x)
x = x[grep("Bmat",x$par),]
x$par = substr(x$par,1,9)

x = dplyr::left_join(x, grid)
x$scenario = paste0(x$obs_sd, ", ", x$pro_sd)

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

g1 = ggplot(x, aes(scenario,mean,group=scenario, col=scenario, fill=scenario)) + 
  geom_boxplot(alpha=0.3,outlier.shape = NA) + 
  geom_point(alpha=0.1)+
  ylab("Estimated B parameter") + 
  xlab("Scenario")+ 
  facet_wrap(~par,scale="free_y") + 
  geom_hline(aes(yintercept = true),col="red",alpha=0.3) + 
  scale_fill_viridis(discrete = TRUE, end=0.7) + 
  scale_color_viridis(discrete = TRUE, end=0.7) + 
  theme_bw() + 
  theme(legend.position='none',strip.background = element_rect(color="black",fill="white")) +
  scale_x_discrete(labels = c('0.2, 0.4' = expression(atop(paste(sigma[obs],"= 0.2"), paste(sigma[pro],"= 0.4"))),
    '0.4, 0.2'   = expression(atop(paste(sigma[obs],"= 0.4"),paste(sigma[pro],"= 0.4")))))
pdf("plots/Figure_3.pdf")
g1
dev.off()

jpeg("plots/Figure_3.jpeg")
g1
dev.off()
