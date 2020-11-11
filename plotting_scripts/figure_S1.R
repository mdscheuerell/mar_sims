library(ggplot)
library(viridis)
fit = readRDS("fit.rds")

fit = as.data.frame(fit)
fitdf = fit[,c("Bmat[1,1]","Bmat[1,2]","SD_proc","SD_obs","lp__")]
names(fitdf) = c("B11","B12","sd_proc","sd_obs","lp")

g1 = ggplot(fitdf, aes(x=sd_obs, y=sd_proc, col=lp)) + 
  geom_hex(bins=50) + scale_fill_viridis_c() + 
  xlab(expression(paste("Obs. error ",sigma))) + 
  ylab(expression(paste("Pro. error ",sigma))) + 
  labs(fill = "Frequency")

g2 = ggplot(fitdf, aes(sd_obs, sd_proc, z=lp)) + 
  stat_summary_hex(aes(z = lp), bins = 50, fun = mean) + scale_fill_viridis_c() + 
  xlab(expression(paste("Obs. error ",sigma))) + 
  ylab(expression(paste("Pro. error ",sigma))) + 
  labs(fill = "Log posterior")

pdf("plots/Figure_S1_variance_tradeoff.pdf")
gridExtra::grid.arrange(g1,g2,nrow=1)
dev.off()

jpeg("plots/Figure_S1_variance_tradeoff.jpeg")
gridExtra::grid.arrange(g1,g2,nrow=1)
dev.off()
