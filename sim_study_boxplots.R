rm(list=ls())
load('results/sim_competing_results_train.RData')
comp <- est_matrix
load('results/sim_SPCDE_results_train.RData')
spcde <- est_matrix
namae <- c('tmax-wasdist','tmax-quantile','tmax-spatcorr',
           'prcp-wasdist','prcp-quantile','prcp-spatcorr',
           'prcp-propzero','crosscorr')

pdf(file = 'sim_results_train.pdf',12,6)
par(mfrow=c(2,4),mgp=c(2.25,0.75,0),mar=c(3,3,1,1))
for(i in c(1,2,3,8,4,5,6,7)){
    metric <- rbind(comp[i,1:2,],spcde[i,2:3,])
    boxplot(t(metric),main=namae[i],cex=0.75,names = c('QM','CCE','SPECD1','SPECD2'), 
            outline=F,col='white')
}
par(mfrow=c(1,1))
dev.off()
