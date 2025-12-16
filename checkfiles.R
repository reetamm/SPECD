folder <- 'fits/sim_SPCDE'
for(i in 1:100){
    folder1 <- paste(folder,i,sep='/')
    allfiles <- list.files(folder1)
    if(length(allfiles)<50)
        print(paste('Folder',i,'has',length(allfiles),'files'))
}
