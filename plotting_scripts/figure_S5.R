library(dplyr)
library(ggplot2)
library(viridis)

grid = readRDS("grid.rds")
grid$iter = seq(1,nrow(grid))

x = readRDS("results/posterior_summaries_C22_1.rds")
x$obsSD = 0.2
x$proSD = 0.2
x$scenario = "0.2, 0.2"
x$par = rownames(x)
x = x[grep("Bmat",x$par),]
x$par = substr(x$par,1,9)
x = dplyr::left_join(x, grid)

x2 = readRDS("results/posterior_summaries_C44_1.rds")
x2$obsSD = 0.4
x2$proSD = 0.4
x2$scenario = "0.4, 0.4"
x2$par = rownames(x2)
x2 = x2[grep("Bmat",x2$par),]
x2$par = substr(x2$par,1,9)
x2 = dplyr::left_join(x2, grid)

x = rbind(x,x2)

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

x = dplyr::filter(x, b_CV==1)

g1 = ggplot(x, aes(scenario,mean,group=scenario, col=scenario, fill=scenario)) + 
  geom_boxplot(alpha=0.4,outlier.shape = NA) + 
  ylab("Estimated B parameter") + 
  xlab("Scenario")+ 
  facet_wrap(~par,scale="free_y") + 
  geom_hline(aes(yintercept = true),col="red",alpha=0.3) + 
  scale_fill_viridis(discrete = TRUE, end=0.7) + 
  scale_color_viridis(discrete = TRUE, end=0.7) + 
  theme_bw() + 
  theme(legend.position='none',strip.background = element_rect(color="black",fill="white")) +
  scale_x_discrete(labels = c('0.2, 0.2' = expression(atop(paste(sigma[obs],"= 0.2"), paste(sigma[pro],"= 0.2"))),
    '0.4, 0.4'   = expression(atop(paste(sigma[obs],"= 0.4"),paste(sigma[pro],"= 0.4")))))
pdf("plots/Figure_S5.pdf")
g1
dev.off()

jpeg("plots/Figure_S5.jpeg")
g1
dev.off()
