
tissue='Brain_Frontal_Cortex_BA9'
tissue='Whole_Blood'  

ut<-read.table(paste0('~/../Dropbox/DansPaper/data/replication/ut_modified_',tissue,'.txt'),header = T,stringsAsFactors = F)
ut<-ut[,-4]
colnames(ut)<-c('geneid','r_ut','p_ut')

xt<-read.table(paste0('~/../Dropbox/DansPaper/data/replication/xt_',tissue,'.txt'),header = T,stringsAsFactors = F)
xt<-xt[,-4]
colnames(xt)<-c('geneid','r_xt','p_xt')

st<-read.table(paste0('~/../Dropbox/DansPaper/data/replication/st_',tissue,'.txt'),header = T,stringsAsFactors = F)
st<-st[,-4]
colnames(st)<-c('geneid','r_st','p_st')

df<-merge(ut,xt,by=1,all=T)
df<-merge(df,st,by=1,all=T)

df[is.na(df[,2]),2]<-0
df[is.na(df[,3]),3]<-1
df[is.na(df[,4]),4]<-0
df[is.na(df[,5]),5]<-1
df[is.na(df[,6]),6]<-0
df[is.na(df[,7]),7]<-1
df[,'r_st']<-ifelse(df[,'r_st']<0,0,df[,'r_st'])
df[,'r_ut']<-ifelse(df[,'r_ut']<0,0,df[,'r_ut'])
df[,'r_xt']<-ifelse(df[,'r_xt']<0,0,df[,'r_xt'])


df$r2_st<-df$r_st^2
df$r2_xt<-df$r_xt^2
df$r2_ut<-df$r_ut^2

df$iG_r_ut<-ifelse((df$r_ut>0.1 & df$p_ut<0.05),1,0)
df$iG_r_xt<-ifelse((df$r_xt>0.1 & df$p_xt<0.05),1,0)
df$iG_r_st<-ifelse((df$r_st>0.1 & df$p_st<0.05),1,0)

df$iG_r2_ut<-ifelse((df$r2_ut>0.01 & df$p_ut<0.05),1,0)
df$iG_r2_xt<-ifelse((df$r2_xt>0.01 & df$p_xt<0.05),1,0)
df$iG_r2_st<-ifelse((df$r2_st>0.01 & df$p_st<0.05),1,0)

sum(df$iG_r_ut)
sum(df$iG_r_xt)
sum(df$iG_r_st)

sum(df$iG_r2_ut)
sum(df$iG_r2_xt)
sum(df$iG_r2_st)


ut_fun<-function(a,b){  #a=st
  if (a==0 & b==0){
    col='grey'
  }else if(a==1 & b==1){
    col='black'
  }else if(a==1 & b==0){
    col=brewer.pal(12, "Paired")[8] #'firebrick2'
  }else if(a==0 & b==1){
    col=brewer.pal(12, "Paired")[2] #'dodgerblue3'
  }
  return(col)
}

xt_fun<-function(a,b){  #a=st
  if (a==0 & b==0){
    col='grey'
  }else if(a==1 & b==1){
    col='black'
  }else if(a==1 & b==0){
    col=brewer.pal(12, "Paired")[8]
  }else if(a==0 & b==1){
    col=brewer.pal(12, "Paired")[4]
  }
}


iG_r_st<- df[df$iG_r_st==1,1]
iG_r_ut<- df[df$iG_r_ut==1,1]
iG_r_xt<- df[df$iG_r_xt==1,1]

iG_r2_st<- df[df$iG_r2_st==1,1]
iG_r2_ut<- df[df$iG_r2_ut==1,1]
iG_r2_xt<- df[df$iG_r2_xt==1,1]



#plot
library (VennDiagram)
library(RColorBrewer)
library("vioplot")

png(paste0('~/../Dropbox/DansPaper/figures/tmp/rep_',tissue,'.png'),width = 2400,height = 600,res=330)
par(mfcol=c(1,4))

#barplot
par(mai=c(1,1,0.5,0.5))
par(mgp=c(3,0.8,0))
par(mar=c(5,5,1,0))
barplot(c(sum(df$iG_r_st),sum(df$iG_r_ut),sum(df$iG_r_xt)), 
        ylab="Number of iGenes",
        ylim = c(0,max(c(sum(df$iG_r_st),sum(df$iG_r_ut),sum(df$iG_r_xt)+800))),
        beside=F,horiz = F,las=2,cex.axis=1,
        col = brewer.pal(12, "Paired")[c(8,2,4)],
        names.arg=c("PrediXcan","UTMOST","XT-SCAN")
)

#venn plot
mgp=c(0,0,0)
plot(c(1),c(1),type="n",xaxt="n",yaxt="n",xlab = '',ylab = '',bty='n')
#venn(x= list('PrediXcan' = iG_r_st,'UTMOST' = iG_r_ut,'XT-SCAN' = iG_r_xt), col="transparent",fill=brewer.pal(12, "Paired")[c(7,1,3)], lwd=0.6, cex=1.5, cat.cex=0.8)


#scatter plot
df$col_ut<-mapply(ut_fun,df$iG_r_st,df$iG_r_ut)
df$col_xt<-mapply(xt_fun,df$iG_r_st,df$iG_r_xt)

#ut v.s. st
par(mar=c(2.5,2.5,1,1),mgp=c(1.5,0.5,0))

plot(-1,-1,xlab = 'PrediXcan',ylab='UTMOST',bty='l',xlim = c(0,1),ylim = c(0,1))
points(df$r_st,df$r_ut,col=df$col_ut,pch=20,cex=0.5)
segments(-1,-1,1,1,col = 'red',lty=2)

#xt v.s. st

plot(-1,-1,xlab = 'PrediXcan',ylab='XT-SCAN',bty='l',xlim = c(0,1),ylim = c(0,1))
points(df$r_st,df$r_xt,col=df$col_xt,pch=20,cex=0.5)
segments(-1,-1,1,1,col = 'red',lty=2)

dev.off()


venn.diagram(x= list('PrediXcan' = iG_r_st,'UTMOST' = iG_r_ut,'XT-SCAN' = iG_r_xt), filename = paste0('~/../Dropbox/DansPaper/figures/tmp/rep_',tissue,'_venn.png'), height = 800, width = 800,,resolution =250, imagetype="png", col="transparent",fill=brewer.pal(12, "Paired")[c(7,1,3)], lwd=0.6, cex=1.5, cat.cex=0.8)






#vio plot


png(paste0('~/../Dropbox/DansPaper/figures/tmp/rep_vio_',tissue,'.png'),width = 500,height = 800,res=150)
vioplot(df$r_st,df$r_ut,df$r_xt,names = c('PrediXcan','UTMOST','XT-SCAN'),col = brewer.pal(12, "Paired")[c(8,2,4)],cex.axis = 1,las=2,ylim=c(0,1),ylab='Pred and Obs correlation r',bty='l',main=tissue)
dev.off()









