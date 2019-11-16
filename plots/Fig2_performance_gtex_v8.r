
#-----performance-----v8----cv-----

tpm<-read.table('~/../Dropbox/data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.49',header = T,stringsAsFactors = F)
tpm[,1]<-sapply(tpm[,1], function(x) strsplit(x,"[.]")[[1]][1])

main_path=paste0('~/../Dropbox/DansPaper/data/performance/v8/')
tissue_list<-dir(main_path)

output<-as.data.frame(matrix(data=NA,ncol=14))
colnames(output)<-c('tissue','n_xt','n_ut','n_st','delta_r2_xt','delta_r2_ut','n_increase_xt','n_increase_ut','delta_r2_xt_p','delta_r2_ut_p','r2_st_median','r2_xt_median','r2_st_mean','r2_xt_mean')


i=1
for (i in 1:length(tissue_list)){
  print(tissue_list[i])
  box<-read.table(paste0(main_path,tissue_list[i]),header = T,stringsAsFactors = F)
  box[is.na(box)] <- 0
  box<-box[which(box[,1] %in% tpm[tpm[,i+2]>0,1]),] #median TPM >0.1
  #print(colnames(tpm[i+2]))
  
  box$delta_xt<-box$r_xt^2-box$r_st^2
  box$delta_ut<-box$r_ut^2-box$r_st^2
  box$delta_xt_p<-(box$r_xt^2-box$r_st^2)/(box$r_st^2+0.001)
  box$delta_ut_p<-(box$r_ut^2-box$r_st^2)/(box$r_st^2+0.001)
  
  output[i,1]<-tissue_list[i]
  
  output[i,2]<-length(which(box$r_xt>0.1 & box$p_xt<0.05))
  output[i,3]<-length(which(box$r_ut>0.1 & box$p_ut<0.05))
  output[i,4]<-length(which(box$r_st>0.1 & box$p_st<0.05))
  output[i,5]<-mean(box$delta_xt)
  output[i,6]<-mean(box$delta_ut)
  output[i,9]<-median(box$delta_xt_p)
  output[i,10]<-median(box$delta_ut_p)
  output[i,11]<-median(box$r_st^2)
  output[i,12]<-median(box$r_xt^2)
  output[i,13]<-mean(box$r_st^2)
  output[i,14]<-mean(box$r_xt^2)
}

output[,7]<-(output[,2]-output[,4])/output[,4]
output[,8]<-(output[,3]-output[,4])/output[,4]

write.table(output,paste0('~/../Dropbox/DansPaper/data/performance/v8_cv.txt'),quote = F,row.names = F,sep='\t')

#-----merge with info-----
info<-read.table('~/../Dropbox/DansPaper/data/info/gtex_v8_info.txt',header = T,stringsAsFactors = F)
info<-info[,-4]
v8<-read.table(paste0('~/../Dropbox/DansPaper/data/performance/v8_cv.txt'),header = T,stringsAsFactors = F)
v8<-merge(info,v8,by=1)


#-----capture-----
tpm<-read.table('~/../Dropbox/data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.49',header = T,stringsAsFactors = F)
tpm[,1]<-sapply(tpm[,1], function(x) strsplit(x,"[.]")[[1]][1])


main_path=paste0('~/../Dropbox/DansPaper/data/performance/v8/')
tissue_list<-dir(main_path)
capture<-as.data.frame(matrix(data=NA,ncol=2))
colnames(capture)<-c('tissue','capture')

i=1
for (i in 1:length(tissue_list)){
  print(tissue_list[i])
  box<-read.table(paste0(main_path,tissue_list[i]),header = T,stringsAsFactors = F)
  box[is.na(box)] <- 0
  box<-box[which(box[,1] %in% tpm[tpm[,i+2]>0,1]),]
  #print(colnames(tpm[i+2]))
  
  box$xt<-ifelse((box$r_xt>0.1 & box$p_xt<0.05),1,0)
  box$st<-ifelse((box$r_st>0.1 & box$p_st<0.05),1,0)
  
  capture[i,1]<-tissue_list[i]
  capture[i,2]<-sum(box[box$st==1,8])/length(which(box$st==1))
}

v8<-merge(v8,capture,by=1)

v8<-v8[order(v8$n,decreasing = T),]
#v8<-v8[order(v8$tissue,decreasing = T),]

write.table(v8,paste0('~/../Dropbox/DansPaper/data/performance/v8_info.txt'),quote = F,row.names = F,sep='\t')

#-----------------------------------------------------------------

png('~/../Dropbox/DansPaper/figures/tmp/v8_performance.png',height = 1200,width = 2000,res=150)
par(mfcol=c(1,3))

#---N_iGenes---
par(mai=c(1,0.1,0.1,0.3))
par(mgp=c(2.5,1,0))
plot(-100,-100,xlim=c(0,max(v8$n_xt)*1.1),ylim=c(0,49),yaxt="n",las=1,xlab='Number of iGenes',ylab=' ',cex.lab=1.5,cex.axis=1.5, bty='l')

#legend
points(200,8,pch=21,bg='white',col='black',cex=2.5)
points(200,6,pch=23,bg='white',col='black',cex=2.5)

text(300,8,'PrediXcan',pos=4,cex=1.2)
text(300,6,'XT-SCAN',pos=4,cex=1.2)

for (i in 1:49){
  segments(v8[i,6],i,v8[i,4],i,col=paste0('#',v8[i,3]))
  points(v8[i,6],i,pch=21,bg=paste0('#',v8[i,3]),cex=2.5)
  points(v8[i,4],i,pch=23,bg=paste0('#',v8[i,3]),cex=2.5)
  
  #mark
  points(18000,i,pch=22,cex=2.5,bg=paste0('#',v8[i,3]))
  
}


#---capture---
par(mai=c(1,2,0.1,0.3))
par(mgp=c(2.5,1,0))

for (i in 1:10){
  v8$tissue<-sub('_',' ',v8$tissue)
}

plot(-100,-100,xlim=c(0,1),ylim=c(0,49),yaxt="n",las=1,xlab='% of PrediXcan iGenes captured by XT-SCAN',ylab=' ',cex.lab=1,cex.axis=1.5,bty='l',xaxt='n')

axis(1,at=seq(0,1,0.2),label=paste0(seq(0,1,0.2)*100,'%'),las=1,cex.axis=1.1) 
axis(2,at=seq(1,49),label=v8$tissue,las=1,cex.axis=1.1) 

for (i in 1:49){
  points(v8[i,17],i,pch=24,bg=paste0('#',v8[i,3]),cex=2.5)
  #mark
  points(0,i,pch=22,bg=paste0('#',v8[i,3]),cex=2.5)
  
  #sample size
  text(0.1,i,v8$n[i],cex = 1)
}

text(0.20,-1.5,'N of samples',pos = 3,cex = 1.2)

#---delta_ipiG---

par(mai=c(1,0.8,0.1,0.3))
par(mgp=c(2.5,1,0))

plot(-100,-100,xlim=c(0,max(v8$n)*1.1),ylim=c(0,max(v8$n_increase_xt)*1.1),yaxt='n',las=1,xlab='Sample size for each tissue',ylab=' ',cex.lab=1.5,cex.axis=1.5,bty='l')
title(ylab = 'Increase in the propotion of iGenes (ΔpiG)', mgp = c(5, 1, 0),cex.lab=1.5)

axis(2,at=seq(0,3,0.1),label=paste0(seq(0,3,0.1)*100,'%'),las=1,cex.axis=1.5) 

points(v8$n,v8$n_increase_xt,pch=21,bg=paste0('#',v8[,3]),cex=2.5)


dev.off()

















#--------r2 compare each tissue-------

v8<-v8[order(v8$tissue),]

for (i in 1:10){
  v8$tissue<-sub('_',' ',v8$tissue)
}

main_path=paste0('~/../Dropbox/DansPaper/data/performance/v8/')
tissue_list<-dir(main_path)

library(vioplot)

png('~/../Dropbox/DansPaper/figures/tmp/v8_performance_each_tissue.png',height = 3000,width = 3000,res=200)

par(mfcol=c(7,7))
par(mgp=c(2.2,1.2,0),mar=c(2.5,3.5,1.5,1))

i=1
for (i in 1:length(tissue_list)){
  print(tissue_list[i])
  box<-read.table(paste0(main_path,tissue_list[i]),header = T,stringsAsFactors = F)
  box[is.na(box)] <- 0
  box<-box[which(box[,1] %in% tpm[tpm[,i+2]>0.1,1]),]
  
  vioplot(box$r_st,box$r_xt,col = paste0('#',v8[i,'rgb']),names=c(paste0('PrediXcan \n',round(median(box$r_st),2)),paste0('XT-SCAN \n',round(median(box$r_xt),2))),main=v8[i,'tissue'],ylab = 'Correlation r')
  
}

dev.off()























