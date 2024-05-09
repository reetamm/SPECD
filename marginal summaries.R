rm(list = ls())
y1.cors.0 = NA
y2.cors.0 = NA
y1y2.cors.0 = NA
y1.cors.1 = NA
y2.cors.1 = NA
y1y2.cors.1 = NA
daysinmonth = c(31,28,31,30,31,30,31,31,30,31,30,31)

cal.data = vector('list',12)
for(mnth in 1:12){
    cal.array = array(dim = c(daysinmonth[mnth]*64,6,25))
    for(loc in 1:25){
        print(paste(mnth,loc))
        envname = paste0('fits/fits_t',mnth,'_l',loc,'.RData')
        load(envname)
        # rm(gcm.long,obs.long)
        # save.image(file = envname)
        y1.cors.0 = c(y1.cors.0,cor(qf.y1.mle.ts[y0==0][-1],qf.y1.mle.ts[y0==0][-n0]))
        y1.cors.1 = c(y1.cors.1,cor(cbind(y11,x11))[1,2])
        y2.cors.0 = c(y2.cors.0,cor(qf.y2.mle.ts[y0==0][-1],qf.y2.mle.ts[y0==0][-n0]))
        y2.cors.1 = c(y2.cors.1,cor(cbind(y21,x21))[1,2])
        y1y2.cors.1 = c(y1y2.cors.1,cor(y1[y0==1],y2[y0==1]))
        y1y2.cors.0 = c(y1y2.cors.0,cor(qf.y1.mle.ts[y0==1],qf.y2.mle.ts[y0==1]))

        cal.array[,1,loc] = y1[y0==0]
        cal.array[,3,loc] = y1[y0==1]
        cal.array[,2,loc] = qf.y1.mle.ts[y0==0]
        cal.array[,4,loc] = y2[y0==0]
        cal.array[,6,loc] = y2[y0==1]
        cal.array[,5,loc] = qf.y2.mle.ts[y0==0]
    }
    cal.data[[mnth]] = cal.array
}

save(y1.cors.0,y1.cors.1,y2.cors.0,y2.cors.1,y1y2.cors.0,y1y2.cors.1,cal.data,
           file = 'time_summary.RData')
maxy = 1.01*max(cal.array[,1:2,mnth,loc])
miny = 0.99*min(cal.array[,1:2,mnth,loc])
plot(1:28,cal.array[,1,mnth,loc],type = 'b',col=1,pch=20,ylim = c(miny,maxy))
lines(1:28,cal.array[,2,mnth,loc],type = 'b',col=2,pch=20)

maxy = 1.01*max(cal.array[,3:4,mnth,loc])
miny = 0.99*min(cal.array[,3:4,mnth,loc])
plot(1:28,cal.array[,3,mnth,loc],type = 'b',col=1,pch=20,ylim = c(miny,maxy))
lines(1:28,cal.array[,4,mnth,loc],type = 'b',col=2,pch=20)
plot(y1y2.cors.1,y1y2.cors.0)
abline(0,1)
head(y1y2.cors.1)
cor(cal.array[,1,1,])

heatmap(cor(cal.array[,3,]))
