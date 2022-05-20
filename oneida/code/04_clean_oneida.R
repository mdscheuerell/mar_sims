library(dplyr)
library(tidyr)

# sources: 
# https://ecommons.cornell.edu/bitstream/handle/1813/67639/Schaffner_cornell_0058O_10463.pdf?sequence=1&isAllowed=y
# https://ecommons.cornell.edu/handle/1813/11228
# https://search.dataone.org/data

# notes: 
# 1. group all daphnia by genus
# 2. include blue greens 
# 3. B[1,4] and B[2,4] ON / B[4,3] OFF
# 4. plot symbols for each station on state trajectory

#Year = 1985
# load data
zoop <- read.csv("oneida/data/cbfs.165.5.csv", stringsAsFactors = FALSE)
zoop = zoop[which(zoop$Year==Year),] # 2015 per Schaffner paper
zoop = zoop[grep("Daphnia", zoop$TAXONOMICGROUP),]
phyto = read.csv("oneida/data/cbfs.32.11.csv") %>%
  dplyr::filter(Year == zoop$Year[1])

phyto_spp = read.csv("oneida/data/cbfs.34.14.data") %>% 
  dplyr::rename(X75_95_TaxaID = TaxonomicCode) %>%
  dplyr::select(X75_95_TaxaID, Division)
phyto_spp$Group = NA
phyto_spp$Group[which(phyto_spp$Division=="Bacillariophyta")] = "Diatoms"
phyto_spp$Group[which(phyto_spp$Division%in%c("Chlorophyta","Chrysophyta"))] = "ChlorChrys"
phyto_spp$Group[which(phyto_spp$Division=="Cyanophyta")] = "BlueGreen"
phyto_spp = dplyr::filter(phyto_spp, !is.na(Group))

phyto = dplyr::left_join(phyto, phyto_spp) %>% 
  dplyr::filter(!is.na(Group))
names(phyto) = tolower(names(phyto))
names(zoop) = tolower(names(zoop))
zoop$taxa = "Daphnia"
phyto = dplyr::rename(phyto, taxa = group) %>% 
  dplyr::select(samplingweek, taxa, density, station)
zoop = dplyr::select(zoop, samplingweek, taxa, density, station)
dat <- rbind(zoop, phyto)

# roll up by sampling week -- across stations
agg = 
  dplyr::arrange(dat, samplingweek, taxa,station) %>%  
  dplyr::filter(station=="Shackelton Point") %>%
  dplyr::group_by(taxa, samplingweek) %>%  
  dplyr::summarise(log_density = log(sum(density,na.rm=T))) %>%
  dplyr::filter(samplingweek %in% 18:46) %>% 
  pivot_wider(names_from=c("samplingweek"),values_from = log_density)
# reorder
idx = sort(as.numeric(colnames(agg)[-1]), index.return=T)
agg = agg[,c(1,idx$ix+1)]
saveRDS(agg,"oneida/data/cleaned_data_agg.rds")


# roll up by sampling week -- by station
dat = dplyr::group_by(dat, taxa, samplingweek,station) %>% 
  dplyr::summarise(log_density = log(sum(density,na.rm=T)))
# drop winter
dat_station = dplyr::filter(dat, samplingweek %in% 18:46,
                    station %in% c("Buoy 109","Buoy 125","Shackelton Point"))
dat_wide = pivot_wider(dat_station,names_from=c("samplingweek"),values_from = log_density)

# reorder
idx = sort(as.numeric(colnames(dat_wide)[-c(1:2)]), index.return=T)
dat_wide = dat_wide[,c(1:2,idx$ix+2)]

saveRDS(dat_wide,"oneida/data/cleaned_data_station.rds")

