library(dplyr)
library(ggplot2)
library(rstan)

fit1 <- readRDS("oneida/output/fitted_model_1site.rds")
fit3 <- readRDS("oneida/output/fitted_model_3site.rds")

pars1 <- rstan::extract(fit1)
pars3 <- rstan::extract(fit3)

# extract Bmat
d <- pars1
names = c("Diatom","Green","Daphnia","Bluegreen")

for(i in 1:dim(d$Bmat)[2]) {
  for(j in 1:dim(d$Bmat)[3]) {
    df <- data.frame("par"=paste0(names[j]," on ",names[i]), "value"=d$Bmat[,i,j])
    if(i+j==2) {
      alldf <- df
    } else {
      alldf <- rbind(alldf,df)
    }
  }
}
alldf$Stations = 1
dat1 <- alldf

d <- pars3
for(i in 1:dim(d$Bmat)[2]) {
  for(j in 1:dim(d$Bmat)[3]) {
    df <- data.frame("par"=paste0(names[j]," on ",names[i]), "value"=d$Bmat[,i,j])
    if(i+j==2) {
      alldf <- df
    } else {
      alldf <- rbind(alldf,df)
    }
  }
}
alldf$Stations = 3
dat <- rbind(dat1, alldf)

dat$par = factor(dat$par, levels = c("Diatom on Diatom", "Green on Diatom", "Bluegreen on Diatom","Daphnia on Diatom",
                                     "Diatom on Green", "Green on Green", "Bluegreen on Green","Daphnia on Green",
                                     "Diatom on Bluegreen", "Green on Bluegreen", "Bluegreen on Bluegreen","Daphnia on Bluegreen",
                                     "Diatom on Daphnia", "Green on Daphnia", "Bluegreen on Daphnia","Daphnia on Daphnia"))

dat$Stations = as.factor(dat$Stations)
g1 <- ggplot(dat, aes(value, col=Stations,fill=Stations)) + 
  geom_density(alpha=0.3,col=NA) + 
  facet_wrap(~par,scale="free") + 
  xlab("") + 
  ylab("Density") + theme_bw() + 
  theme(strip.background =element_rect(fill="white")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.x = element_text(size = 7))
pdf("oneida/plots/Oneida_Bmatrix.pdf")
g1
dev.off()

write.csv(apply(pars1$Bmat,c(2,3),sd),"Bmatrix_sd_1site.csv")
write.csv(apply(pars3$Bmat,c(2,3),sd),"Bmatrix_sd_3site.csv")
##############################################################
# make plots of states
xx1 <- broom.mixed::tidy(fit1,"xx",conf.int=TRUE)
xx3 <- broom.mixed::tidy(fit3,"xx",conf.int=TRUE)

xx1$species = rep(1:4, dim(xx1)[1]/4)
xx3$species = rep(1:4, dim(xx3)[1]/4)
xx1$time = sort(rep(1:(dim(xx1)[1]/4),4)) + as.numeric(colnames(dat_wide)[3]) - 1
xx3$time = sort(rep(1:(dim(xx3)[1]/4),4)) + as.numeric(colnames(dat_wide)[3]) - 1
xx1$Stations = 1
xx3$Stations = 3
xx <- rbind(xx1,xx3)
xx$Stations = as.factor(xx$Stations)
xx$Species = c("Cyanophyta","Chlorophyta & Chrysophyta","Daphnia","Bacillariophyta")[xx$species]

station_dat1 = data.frame(time = sort(rep(as.numeric(colnames(dat_station1)),4)),
                            Species = rep(c("Cyanophyta","Chlorophyta & Chrysophyta","Daphnia","Bacillariophyta"),ncol(dat_station1)),
                          Stations = "Buoy 109",
                          estimate=c(dat_station1))
station_dat2 = data.frame(time = sort(rep(as.numeric(colnames(dat_station2)),4)),
                          Species = rep(c("Cyanophyta","Chlorophyta & Chrysophyta","Daphnia","Bacillariophyta"),ncol(dat_station2)),
                          Stations = "Buoy 125",
                          estimate=c(dat_station2))
station_dat3 = data.frame(time = sort(rep(as.numeric(colnames(dat_station3)),4)),
                          Species = rep(c("Cyanophyta","Chlorophyta & Chrysophyta","Daphnia","Bacillariophyta"),ncol(dat_station3)),
                          Stations = "Shackelton Point",
                          estimate=c(dat_station3))
station_dat1 = dplyr::filter(station_dat1, time<43)
station_dat2 = dplyr::filter(station_dat2, time<43)
station_dat3 = dplyr::filter(station_dat3, time<43)

station_dat1$Name = "Buoy 109"
station_dat2$Name = "Buoy 125"
station_dat3$Name = "Shackelton Point"
station_dat = rbind(station_dat1, station_dat2, station_dat3)
xx$Name = NA

xx$Species = factor(xx$Species, levels = c("Bacillariophyta","Chlorophyta & Chrysophyta","Cyanophyta","Daphnia"))
g2 <- ggplot(xx, aes(time, estimate, col=Stations, fill=Stations)) + 
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.5, col=NA) + 
  facet_wrap(~Species,scale="free",ncol=1) + 
  geom_line(size=0.3) + 
  ylab("Estimated log density") + 
  xlab("") + 
  theme_bw() + 
  #coord_cartesian(xlim=c(19,41)) + 
  theme(strip.background =element_rect(fill="white")) +
  geom_point(data = station_dat,aes(fill=NA,time,estimate,shape=Name),col="black") + 
  scale_shape_manual(values=c(3, 4, 6))
pdf("Oneida_state_estimates.pdf")
g2
dev.off()
  #geom_point(data = station_dat1, aes(time,estimate),col="black",fill=NA,shape=3) +
  #geom_point(data = station_dat2, aes(time,estimate),col="black",fill=NA,shape=4) + 
  #geom_point(data = station_dat3, aes(time,estimate),col="black",fill=NA,shape=6)
  #geom_point(data = station_dat, aes(time,estimate))

g3 <- ggplot(xx, aes(time, std.error, col=Stations)) + 
  geom_line() + 
  #geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=0.5) + 
  facet_wrap(~Species,scale="free",ncol=1) + 
  ylab("Standard error") + 
  xlab("") + 
  theme_bw() + 
  theme(strip.background =element_rect(fill="white"))
pdf("Fig_S9.pdf")
g3
dev.off()
library(ggpubr)
library(grid)
g <- ggarrange(g2,g3, ncol=2, common.legend = TRUE, legend="top")
g <- annotate_figure(g,
                bottom = textGrob("Sampling week", vjust = -1, gp = gpar(cex = 1))) + 
  ggtitle(paste(Year))


pdf("oneida/plots/CLBS_results.pdf")
g
g1
dev.off()