
options(stringsAsFactors = F)

library(RSQLite)
library(LSD)

#---------------------------------------


#---UTMOST---original---Whole_Blood

tissue='Whole_Blood'
model='ut_original'

#gtex
con <- dbConnect(RSQLite::SQLite(), dbname=paste0('~/../Dropbox/DansPaper/data/db/',model,'_',tissue,'.db')) #establish connections
gtex = dbReadTable(con,"extra")
dbDisconnect(con) #disconnect

gtex<-gtex[,c(1,4)]
colnames(gtex)<-c('genename','r2_gtex')

#rep
rep<-read.table(paste0('~/../Dropbox/DansPaper/data/replication/',model,'_',tissue,'.txt'),header = T)
rep$r<-ifelse(rep$r<0,0,rep$r)
rep$r2_rep<-rep$r^2
rep<-rep[,c(1,4)]

#merge
df<-merge(gtex,rep,by=1)

#plot
png(paste0('~/../Dropbox/DansPaper/figures/tmp/',tissue,'_',model,'.png'),width = 500,height = 1200,res = 150)
par(mfrow=c(2,1))

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,1),ylim = c(0,1),main='UTMOST Original',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (GEUVADIS)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.1,0.1,0.1,col='forestgreen',lty=2)
segments(0.1,0.1,0.1,-1,col='forestgreen',lty=2)

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,.1),ylim = c(0,.1),main='UTMOST Original (zoom in)',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (GEUVADIS)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.1,0.1,0.1,col='forestgreen',lty=2)
segments(0.1,0.1,0.1,-1,col='forestgreen',lty=2)

dev.off()





#---UTMOST---original---Brain_Frontal_Cortex_BA9

tissue='Brain_Frontal_Cortex_BA9'
model='ut_original'

#gtex
con <- dbConnect(RSQLite::SQLite(), dbname=paste0('~/../Dropbox/DansPaper/data/db/',model,'_',tissue,'.db')) #establish connections
gtex = dbReadTable(con,"extra")
dbDisconnect(con) #disconnect

gtex<-gtex[,c(1,4)]
colnames(gtex)<-c('genename','r2_gtex')

#rep
rep<-read.table(paste0('~/../Dropbox/DansPaper/data/replication/',model,'_',tissue,'.txt'),header = T)
rep$r<-ifelse(rep$r<0,0,rep$r)
rep$r2_rep<-rep$r^2
rep<-rep[,c(1,4)]

#merge
df<-merge(gtex,rep,by=1)

#plot
png(paste0('~/../Dropbox/DansPaper/figures/tmp/',tissue,'_',model,'.png'),width = 500,height = 1200,res = 150)
par(mfrow=c(2,1))

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,1),ylim = c(0,1),main='UTMOST Original',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (PsychENCODE)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.3,0.3,0.3,col='forestgreen',lty=2)
segments(0.3,0.3,0.3,-1,col='forestgreen',lty=2)

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,.3),ylim = c(0,.3),main='UTMOST Original (zoom in)',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (PsychENCODE)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.3,0.3,0.3,col='forestgreen',lty=2)
segments(0.3,0.3,0.3,-1,col='forestgreen',lty=2)

dev.off()











#---UTMOST---modified---Whole_Blood

tissue='Whole_Blood'
model='ut_modified'

#gtex
con <- dbConnect(RSQLite::SQLite(), dbname=paste0('~/../Dropbox/DansPaper/data/db/',model,'_',tissue,'.db')) #establish connections
gtex = dbReadTable(con,"extra")
dbDisconnect(con) #disconnect

gtex<-gtex[,c(1,4)]
colnames(gtex)<-c('genename','r2_gtex')

#rep
rep<-read.table(paste0('~/../Dropbox/DansPaper/data/replication/',model,'_',tissue,'.txt'),header = T)
rep$r<-ifelse(rep$r<0,0,rep$r)
rep$r2_rep<-rep$r^2
rep<-rep[,c('gene','r2_rep')]

#merge
df<-merge(gtex,rep,by=1)

#plot
png(paste0('~/../Dropbox/DansPaper/figures/tmp/',tissue,'_',model,'.png'),width = 500,height = 1200,res = 150)
par(mfrow=c(2,1))

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,1),ylim = c(0,1),main='UTMOST Modified',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (GEUVADIS)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.1,0.1,0.1,col='forestgreen',lty=2)
segments(0.1,0.1,0.1,-1,col='forestgreen',lty=2)

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,.1),ylim = c(0,.1),main='UTMOST Modified (zoom in)',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (GEUVADIS)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.1,0.1,0.1,col='forestgreen',lty=2)
segments(0.1,0.1,0.1,-1,col='forestgreen',lty=2)

dev.off()





#---UTMOST---original---Brain_Frontal_Cortex_BA9

tissue='Brain_Frontal_Cortex_BA9'
model='ut_modified'

#gtex
con <- dbConnect(RSQLite::SQLite(), dbname=paste0('~/../Dropbox/DansPaper/data/db/',model,'_',tissue,'.db')) #establish connections
gtex = dbReadTable(con,"extra")
dbDisconnect(con) #disconnect

gtex<-gtex[,c(1,4)]
colnames(gtex)<-c('genename','r2_gtex')

#rep
rep<-read.table(paste0('~/../Dropbox/DansPaper/data/replication/',model,'_',tissue,'.txt'),header = T)
rep$r<-ifelse(rep$r<0,0,rep$r)
rep$r2_rep<-rep$r^2
rep<-rep[,c('gene','r2_rep')]

#merge
df<-merge(gtex,rep,by=1)

#plot
png(paste0('~/../Dropbox/DansPaper/figures/tmp/',tissue,'_',model,'.png'),width = 500,height = 1200,res = 150)
par(mfrow=c(2,1))

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,1),ylim = c(0,1),main='UTMOST Modified',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (PsychENCODE)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.3,0.3,0.3,col='forestgreen',lty=2)
segments(0.3,0.3,0.3,-1,col='forestgreen',lty=2)

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,.3),ylim = c(0,.3),main='UTMOST Modified (zoom in)',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (PsychENCODE)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.3,0.3,0.3,col='forestgreen',lty=2)
segments(0.3,0.3,0.3,-1,col='forestgreen',lty=2)

dev.off()
















#---PrediXcan---Whole_Blood

tissue='Whole_Blood'
model='st'

#gtex
con <- dbConnect(RSQLite::SQLite(), dbname=paste0('~/../Dropbox/DansPaper/data/db/',model,'_',tissue,'.db')) #establish connections
gtex = dbReadTable(con,"extra")
dbDisconnect(con) #disconnect

gtex<-gtex[,c(1,4)]
colnames(gtex)<-c('genename','r2_gtex')

#rep
rep<-read.table(paste0('~/../Dropbox/DansPaper/data/replication/',model,'_',tissue,'.txt'),header = T)
rep$r<-ifelse(rep$r<0,0,rep$r)
rep$r2_rep<-rep$r^2
rep<-rep[,c('gene','r2_rep')]

#merge
df<-merge(gtex,rep,by=1)

#plot
png(paste0('~/../Dropbox/DansPaper/figures/tmp/',tissue,'_',model,'.png'),width = 500,height = 1200,res = 150)
par(mfrow=c(2,1))

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,1),ylim = c(0,1),main='PrediXcan',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (GEUVADIS)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.1,0.1,0.1,col='forestgreen',lty=2)
segments(0.1,0.1,0.1,-1,col='forestgreen',lty=2)

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,.1),ylim = c(0,.1),main='PrediXcan (zoom in)',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (GEUVADIS)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.1,0.1,0.1,col='forestgreen',lty=2)
segments(0.1,0.1,0.1,-1,col='forestgreen',lty=2)

dev.off()





#---PrediXcan---Brain_Frontal_Cortex_BA9

tissue='Brain_Frontal_Cortex_BA9'
model='st'

#gtex
con <- dbConnect(RSQLite::SQLite(), dbname=paste0('~/../Dropbox/DansPaper/data/db/',model,'_',tissue,'.db')) #establish connections
gtex = dbReadTable(con,"extra")
dbDisconnect(con) #disconnect

gtex<-gtex[,c(1,4)]
colnames(gtex)<-c('genename','r2_gtex')

#rep
rep<-read.table(paste0('~/../Dropbox/DansPaper/data/replication/',model,'_',tissue,'.txt'),header = T)
rep$r<-ifelse(rep$r<0,0,rep$r)
rep$r2_rep<-rep$r^2
rep<-rep[,c('gene','r2_rep')]

#merge
df<-merge(gtex,rep,by=1)

#plot
png(paste0('~/../Dropbox/DansPaper/figures/tmp/',tissue,'_',model,'.png'),width = 500,height = 1200,res = 150)
par(mfrow=c(2,1))

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,1),ylim = c(0,1),main='PrediXcan',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (PsychENCODE)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.3,0.3,0.3,col='forestgreen',lty=2)
segments(0.3,0.3,0.3,-1,col='forestgreen',lty=2)

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,.3),ylim = c(0,.3),main='PrediXcan (zoom in)',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (PsychENCODE)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.3,0.3,0.3,col='forestgreen',lty=2)
segments(0.3,0.3,0.3,-1,col='forestgreen',lty=2)

dev.off()

















#---XT-SCAN---Whole_Blood

tissue='Whole_Blood'
model='xt'

#gtex
con <- dbConnect(RSQLite::SQLite(), dbname=paste0('~/../Dropbox/DansPaper/data/db/',model,'_',tissue,'.db')) #establish connections
gtex = dbReadTable(con,"extra")
dbDisconnect(con) #disconnect

gtex<-gtex[,c(1,4)]
colnames(gtex)<-c('genename','r2_gtex')

#rep
rep<-read.table(paste0('~/../Dropbox/DansPaper/data/replication/',model,'_',tissue,'.txt'),header = T)
rep$r<-ifelse(rep$r<0,0,rep$r)
rep$r2_rep<-rep$r^2
rep<-rep[,c('gene','r2_rep')]

#merge
df<-merge(gtex,rep,by=1)

#plot
png(paste0('~/../Dropbox/DansPaper/figures/tmp/',tissue,'_',model,'.png'),width = 500,height = 1200,res = 150)
par(mfrow=c(2,1))

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,1),ylim = c(0,1),main='XT-SCAN',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (GEUVADIS)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.1,0.1,0.1,col='forestgreen',lty=2)
segments(0.1,0.1,0.1,-1,col='forestgreen',lty=2)

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,.1),ylim = c(0,.1),main='XT-SCAN (zoom in)',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (GEUVADIS)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.1,0.1,0.1,col='forestgreen',lty=2)
segments(0.1,0.1,0.1,-1,col='forestgreen',lty=2)

dev.off()





#---XT-SCAN---Brain_Frontal_Cortex_BA9

tissue='Brain_Frontal_Cortex_BA9'
model='xt'

#gtex
con <- dbConnect(RSQLite::SQLite(), dbname=paste0('~/../Dropbox/DansPaper/data/db/',model,'_',tissue,'.db')) #establish connections
gtex = dbReadTable(con,"extra")
dbDisconnect(con) #disconnect

gtex<-gtex[,c(1,4)]
colnames(gtex)<-c('genename','r2_gtex')

#rep
rep<-read.table(paste0('~/../Dropbox/DansPaper/data/replication/',model,'_',tissue,'.txt'),header = T)
rep$r<-ifelse(rep$r<0,0,rep$r)
rep$r2_rep<-rep$r^2
rep<-rep[,c('gene','r2_rep')]

#merge
df<-merge(gtex,rep,by=1)

#plot
png(paste0('~/../Dropbox/DansPaper/figures/tmp/',tissue,'_',model,'.png'),width = 500,height = 1200,res = 150)
par(mfrow=c(2,1))

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,1),ylim = c(0,1),main='XT-SCAN',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (PsychENCODE)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.3,0.3,0.3,col='forestgreen',lty=2)
segments(0.3,0.3,0.3,-1,col='forestgreen',lty=2)

heatscatter(df$r2_gtex,df$r2_rep,pch='.',xlim = c(0,.3),ylim = c(0,.3),main='XT-SCAN (zoom in)',xlab='r2 in the training set (GTEx)',ylab='r2 in the test set (PsychENCODE)',cexplot = 2)
segments(-1,-1,2,2,col='black',lty=2)
segments(-1,0.3,0.3,0.3,col='forestgreen',lty=2)
segments(0.3,0.3,0.3,-1,col='forestgreen',lty=2)

dev.off()










