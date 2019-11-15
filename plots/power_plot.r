#power plot


info<-read.table('~/../Dropbox/DansPaper/data/info/gtex_v8_info.txt',header = T,stringsAsFactors = F)
info$x0.10_st<-NA
info$x0.30_st<-NA
info$x0.50_st<-NA
info$x0.10_ct<-NA
info$x0.30_ct<-NA
info$x0.50_ct<-NA

tissue_list<-info$tissue

model='st'
i=1
for (model in c('st','xt')){
  for (i in 1:length(tissue_list)){
    df<-read.table(paste0('~/../Dropbox/DansPaper/data/simulation/power/',model,'/',tissue_list[i],'.txt'),header = T,stringsAsFactors = F)
    df<-df[!duplicated(df$gene),]
    
    if (model=='st'){
      info[i,4]<-mean(df$X0.10)
      info[i,5]<-mean(df$X0.30)  
      info[i,6]<-mean(df$X0.50)  
    }else{
      info[i,7]<-mean(df$X0.10)
      info[i,8]<-mean(df$X0.30) 
      info[i,9]<-mean(df$X0.50)  
    }
    
  }
}

info<-info[order(info$n*-1),]

for (i in 1:10){
  info[,1]<-sub('_',' ',info[,1])
}

write.table(info,'~/../Dropbox/DansPaper/data/simulation/power/power.txt',quote = F,sep='\t',row.names = F)


png('~/../Dropbox/DansPaper/figures/tmp/power.png',height = 1600,width = 1500,res=150)
par(mfcol=c(1,1))

#N of iGenes
par(mai=c(1,4,0.1,0.3))
par(mgp=c(2.5,1,0))
plot(-100,-100,xlim=c(0,1),ylim=c(-2,49),yaxt="n",las=1,xlab='Power',ylab=' ',cex.lab=1.8,cex.axis=1.8,bty='l')

#marker
segments(0.08,-0.5,0.13,-0.5)
segments(0.08,-1.5,0.13,-1.5)
segments(0.08,-2.5,0.13,-2.5)

points(0.08,-0.5,pch=24,cex=1.4,bg='white')
points(0.13,-0.5,pch=24,cex=1.4,bg='black')
text(0.16,-0.5,'α_g=0.1 (PrediXcan & XT-SCAN)',pos=4,cex=1)

points(0.08,-1.5,pch=23,cex=1.4,bg='white')
points(0.13,-1.5,pch=23,cex=1.4,bg='black')
text(0.16,-1.5,'α_g=0.3 (PrediXcan & XT-SCAN)',pos=4,cex=1)

points(0.08,-2.5,pch=21,cex=1.4,bg='white')
points(0.13,-2.5,pch=21,cex=1.4,bg='black')
text(0.16,-2.5,'α_g=0.5 (PrediXcan & XT-SCAN)',pos=4,cex=1)


axis(2,at=seq(1,49),label=info$tissue,las=1,cex.axis=1.2) 

for (i in 1:49){
  print(i)
  
  segments(info$x0.10_st[i],i+0.3,info$x0.10_ct[i],i+0.3,col=paste0('#',info$rgb[i]))
  segments(info$x0.30_st[i],i,info$x0.30_ct[i],i,col=paste0('#',info$rgb[i]))
  segments(info$x0.50_st[i],i-0.3,info$x0.50_ct[i],i-0.3,col=paste0('#',info$rgb[i]))
  
  points(info$x0.10_st[i],i+0.3,bg='white',pch=24,cex=1.4)
  points(info$x0.10_ct[i],i+0.3,bg=paste0('#',info$rgb[i]),pch=24,cex=1.4)
  
  points(info$x0.30_st[i],i,bg='white',pch=23,cex=1.4)
  points(info$x0.30_ct[i],i,bg=paste0('#',info$rgb[i]),pch=23,cex=1.4)
  
  points(info$x0.50_st[i],i-0.3,bg='white',pch=21,cex=1.4)
  points(info$x0.50_ct[i],i-0.3,bg=paste0('#',info$rgb[i]),pch=21,cex=1.4)
  
  #mark
  points(-0.02,i,pch=22,cex=2.5,bg=paste0('#',info$rgb[i]))
}

dev.off()
