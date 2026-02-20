library(ggplot2)
library(lubridate)

# Obs prcp
obs.long <- read.table('anomaly data/selected_nclim_prcp_anomaly.csv',sep = ',',header = T)
# head(obs.long)
nrow(obs.long)/length(unique(obs.long$lat*obs.long$lon))
obs.months = month(obs.long[,2])
obs.years = year(obs.long[,2])
summary(obs.long$X)
numdays = c(31,28,31,30,31,30,31,31,30,31,30,31)
obs_prcp_1 <- obs.long[obs.months==1 & obs.years==2014,c(3,4,6)]
obs_prcp_1_agg <- aggregate(obs_prcp_1,by=list(obs_prcp_1$lat, obs_prcp_1$lon),FUN=mean)
ggplot(obs_prcp_1_agg,aes(x=lon,y=lat,fill=prcp)) +
    geom_tile() +
    scale_fill_viridis_c(option = 'H') +
    coord_fixed() + theme_bw() + 
    labs(fill='PRCP',x='',y='',title='Jan 2014  nClimGrid PRCP anomalies (mm)') +
    theme(plot.title = element_text(hjust=.5))
ggsave('obs_prcp_1.pdf',width = 4,height = 3)

obs_prcp_7 <- obs.long[obs.months==7 & obs.years==2014,c(3,4,6)]
obs_prcp_7_agg <- aggregate(obs_prcp_7,by=list(obs_prcp_7$lat, obs_prcp_7$lon),FUN=mean)
ggplot(obs_prcp_7_agg,aes(x=lon,y=lat,fill=prcp)) +
    geom_tile() +
    scale_fill_viridis_c(option = 'H') +
    coord_fixed() + theme_bw() + 
    labs(fill='PRCP',x='',y='',title='Jul 2014 nClimGrid PRCP anomalies (mm)') +
    theme(plot.title = element_text(hjust=.5))
ggsave('obs_prcp_7.pdf',width = 4,height = 3)

###### tmax obs

obs.long <- read.table('anomaly data/selected_nclim_tmax_anomaly.csv',sep = ',',header = T)
# head(obs.long)
nrow(obs.long)/length(unique(obs.long$lat*obs.long$lon))
obs.months = month(obs.long[,2])
obs.years = year(obs.long[,2])
summary(obs.long$X)
numdays = c(31,28,31,30,31,30,31,31,30,31,30,31)
obs_tmax_1 <- obs.long[obs.months==1 & obs.years==2014,c(3,4,6)]
obs_tmax_1_agg <- aggregate(obs_tmax_1,by=list(obs_tmax_1$lat, obs_tmax_1$lon),FUN=mean)
ggplot(obs_tmax_1_agg,aes(x=lon,y=lat,fill=tmax)) +
    geom_tile() +
    scale_fill_viridis_c(option = 'H') +
    coord_fixed() + theme_bw() + 
    labs(fill='TMAX',x='',y='',title='Jan 2014  nClimGrid TMAX anomalies (K)') +
    theme(plot.title = element_text(hjust=.5))
ggsave('obs_tmax_1.pdf',width = 4,height = 3)

obs_tmax_7 <- obs.long[obs.months==7 & obs.years==2014,c(3,4,6)]
obs_tmax_7_agg <- aggregate(obs_tmax_7,by=list(obs_tmax_7$lat, obs_tmax_7$lon),FUN=mean)
ggplot(obs_tmax_7_agg,aes(x=lon,y=lat,fill=tmax)) +
    geom_tile() +
    scale_fill_viridis_c(option = 'H') +
    coord_fixed() + theme_bw() + 
    labs(fill='TMAX',x='',y='',title='Jul 2014 nClimGrid TMAX anomalies (K)') +
    theme(plot.title = element_text(hjust=.5))
ggsave('obs_tmax_7.pdf',width = 4,height = 3)
####################
# gcm prcp
region = 'SE'
method = 'MLE'
gcm.long = read.csv(paste0('data/',region,'_gcm_data.csv'))
# obs.long = read.csv(paste0('data/',region,'_obs_data.csv'))
gcm.long$lon <- gcm.long$lon - 360
# head(gcm.long)
nrow(gcm.long)/length(unique(gcm.long$lat*gcm.long$lon))
gcm.months = month(gcm.long[,1])
gcm.years = year(gcm.long[,1])

numdays = c(31,28,31,30,31,30,31,31,30,31,30,31)
gcm_prcp_1 <- gcm.long[gcm.months==1 & gcm.years==2014,2:4]
gcm_prcp_1_agg <- aggregate(gcm_prcp_1,by=list(gcm_prcp_1$lat, gcm_prcp_1$lon),FUN=mean)
ggplot(gcm_prcp_1_agg,aes(x=lon,y=lat,fill=pr)) +
    geom_tile() +
    scale_fill_viridis_c(option = 'H') +
    coord_fixed() + theme_bw() + 
    labs(fill='PRCP',x='',y='',title='Jan 2014 GCM PRCP (mm)') +
    theme(plot.title = element_text(hjust=.5))
ggsave('gcm_prcp_1.pdf',width = 4,height = 3)

gcm_prcp_7 <- gcm.long[gcm.months==7 & gcm.years==2014,2:4]
gcm_prcp_7_agg <- aggregate(gcm_prcp_7,by=list(gcm_prcp_7$lat, gcm_prcp_7$lon),FUN=mean)
ggplot(gcm_prcp_7_agg,aes(x=lon,y=lat,fill=pr)) +
    geom_tile() +
    scale_fill_viridis_c(option = 'H') +
    coord_fixed() + theme_bw() + 
    labs(fill='PRCP',x='',y='',title='Jul 2014 GCM PRCP (mm)') +
    theme(plot.title = element_text(hjust=.5))
ggsave('gcm_prcp_7.pdf',width = 4,height = 3)

###### tmax gcm

region = 'SE'
method = 'MLE'
gcm.long = read.csv(paste0('data/',region,'_gcm_data.csv'))
# obs.long = read.csv(paste0('data/',region,'_obs_data.csv'))
gcm.long$lon <- gcm.long$lon - 360
head(gcm.long)
nrow(gcm.long)/length(unique(gcm.long$lat*gcm.long$lon))
gcm.months = month(gcm.long[,1])
gcm.years = year(gcm.long[,1])
# summary(gcm.long$X)
numdays = c(31,28,31,30,31,30,31,31,30,31,30,31)
gcm_tmax_1 <- gcm.long[gcm.months==1 & gcm.years==2014,c(2,3,5)]
gcm_tmax_1_agg <- aggregate(gcm_tmax_1,by=list(gcm_tmax_1$lat, gcm_tmax_1$lon),FUN=mean)
ggplot(gcm_tmax_1_agg,aes(x=lon,y=lat,fill=tmax)) +
    geom_tile() +
    scale_fill_viridis_c(option = 'H') +
    coord_fixed() + theme_bw() + 
    labs(fill='TMAX',x='',y='',title='Jan 2014 GCM TMAX (K)') +
    theme(plot.title = element_text(hjust=.5))
ggsave('gcm_tmax_1.pdf',width = 4,height = 3)

gcm_tmax_7 <- gcm.long[gcm.months==7 & gcm.years==2014,c(2,3,5)]
gcm_tmax_7_agg <- aggregate(gcm_tmax_7,by=list(gcm_tmax_7$lat, gcm_tmax_7$lon),FUN=mean)
ggplot(gcm_tmax_7_agg,aes(x=lon,y=lat,fill=tmax)) +
    geom_tile() +
    scale_fill_viridis_c(option = 'H') +
    coord_fixed() + theme_bw() + 
    labs(fill='TMAX',x='',y='',title='Jul 2014 GCM TMAX (K)') +
    theme(plot.title = element_text(hjust=.5))
ggsave('gcm_tmax_7.pdf',width = 4,height = 3)

