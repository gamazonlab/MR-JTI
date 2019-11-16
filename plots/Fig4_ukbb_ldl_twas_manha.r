#st xt

library(RColorBrewer)
library(ggrepel)
library(ggplot2)
#display.brewer.all(type = "qual")

#---colors---
dark_col<-brewer.pal(12, "Paired")[c(8,2,4)]  #orange blue green
light_col<-brewer.pal(12, "Paired")[c(7,1,3)]  #orange blue green
names(dark_col)<-names(light_col)<-c('st','ut','xt')

#---load_annotation---
gencode<-read.table('~/../Dropbox/annotation/encode/gencode.v32.GRCh38.txt',header = T,stringsAsFactors = F)
gencode<-gencode[,c('geneid','chr','left','right','genename')]


#---load_associations---
st<-read.csv('~/../Dropbox/DansPaper/data/apply_gwas/LDLq_st_Liver.csv',header = T,stringsAsFactors = F)
xt<-read.csv('~/../Dropbox/DansPaper/data/apply_gwas/LDLq_xt_Liver.csv',header = T,stringsAsFactors = F)

pvalue.extreme <- function(z) {
  log.pvalue <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
  log10.pvalue <- log.pvalue/log(10) ## from natural log to log10
  mantissa <- 10^(log10.pvalue %% 1)
  exponent <- log10.pvalue %/% 1
  
  logp<-log10.pvalue*-1
  return(logp)
  ##return(sprintf("p value is %1.2f times 10^(%d)",mantissa,exponent))
}

st$plog<-pvalue.extreme(st$zscore)
xt$plog<-pvalue.extreme(xt$zscore)


st<-st[,c(1,ncol(st))];colnames(st)<-c('geneid','plog')
st$know_col=dark_col[1];st$add_col=light_col[1];st$order=3;st$model='PrediXcan'
xt<-xt[,c(1,ncol(xt))];colnames(xt)<-c('geneid','plog')
xt$know_col=dark_col[3];xt$add_col=light_col[3];xt$order=1;xt$model='XT-SCAN'

df<-rbind(st,xt)
df<-merge(gencode,df,by='geneid')

df$pos<-(df$left+df$right)/2

df$chr=as.numeric(sapply(df$chr,function(x) sub('^...','',x)))


#---generate_pos---
chr_box<-as.data.frame(matrix(data=NA,nrow=22,ncol=3))
chr_df<-df[df$chr==1,]
chr_box[1,1]<-max(chr_df$pos)
chr_box[1,2]<-0
chr_box[1,3]<-chr_box[1,1]/2

for (i in 2:22){
  chr_df<-df[df$chr==i,]
  chr_box[i,1]<-max(chr_df$pos)
  chr_box[i,2]<-max(chr_df$pos)+chr_box[i-1,2]
  df[which(df$chr==i),'pos']<-df[which(df$chr==i),'pos']+chr_box[i,2]
  chr_box[i,3]<-chr_box[i,2]+chr_box[i,1]/2
}

#---rm_non_sig_genes---
df<-df[df$plog>log(0.05,10)*-1,]


#---load know genes---
know<-read.csv('~/../Dropbox/XTSCAN/glgc/knowgenes.csv',header = T,stringsAsFactors = F)
df[!(df$genename %in% know$genename),'genename']<-'ZZZ'
df$col=ifelse(df$genename %in% know$genename,df$know_col,df$add_col)

#---order---
df<-df[order(df$genename,df$order,decreasing = T),]
df$genename<-ifelse(df$genename=='ZZZ',NA,df$genename)

#---name---
df$group<-ifelse(!(df$genename %in% know$genename),paste0('Additional ',df$model),paste0('Known ',df$model))
df$group<-as.factor(df$group)

#---reorder legend---
df$group <- factor(df$group, levels = rev(levels(df$group)))

#---ggplot---

png('~/../Dropbox/DansPaper/figures/tmp/ukbb_ldl_man.png',width = 2200,height = 1000,res=250)

set.seed(1)
ggplot(df, aes(x=pos, y=plog, label=genename)) + 
  scale_x_continuous(breaks=chr_box$V3, labels = c(as.character(seq(1,15)),' ','17','','19','','21','')) +  #x axis
  
  scale_y_continuous(breaks=c(1,2,3,5,10,20,30,50,100,200,300,500),trans='log10') +  #x axis
  
  #scale_y_continuous(trans='log10') +  #y axis log
  geom_point(data = df, aes(x = pos, y = plog, color = group, shape=group,fill=group)) + #points
  scale_fill_manual(values = c("Known PrediXcan"=as.character(dark_col[1]), 'Known XT-SCAN'=as.character(dark_col[3]),"Additional PrediXcan"=as.character(light_col[1]), 'Additional XT-SCAN'=as.character(light_col[3]))) +
  scale_color_manual(values = c("Known PrediXcan"='black', 'Known XT-SCAN'='black',"Additional PrediXcan"=as.character(light_col[1]),'Additional XT-SCAN'=as.character(light_col[3]))) +
  scale_shape_manual(values = c("Known PrediXcan"=21, 'Known XT-SCAN'=21,"Additional PrediXcan"=20, 'Additional XT-SCAN'=20)) +
  
  labs(x = "Chromosome", y = "-log(P)", title = "") +
  theme_bw() +  #rm background
  theme(panel.grid =element_blank()) +  #rm grids
  theme(legend.position= "right",legend.title = element_blank()) + #set legend
  geom_label_repel( #non overlapped labels
    size=2, 
    fill=rgb(255, 255, 255, 200, maxColorValue=255),
    colour = df$col,  
    direction='y', #only work on Y-axis
    nudge_x=1.5e8, #shift to the right
    segment.alpha = 0.2,  #transparent of segment
    min.segment.length = 0.5,
    fontface='bold.italic'
  )

dev.off()







#bar

png('~/../Dropbox/DansPaper/figures/tmp/ukbb_bar.png',width = 600,height = 1000,res=180)
par(mfcol=c(1,2),mar=c(5,4,3,1),mgp=c(2.5,0.75,0.25))

#known
known<-as.numeric(table(df$group))[c(2,1)]

bar<-barplot(known, border=F, 
             names.arg=c('PrediXcan','XT-SCAN'),cex.names = 0.8,
             las=2, 
             col=c(dark_col['st'],dark_col['xt']), 
             main=paste0("Known \n genes"),cex.main=1,
             ylim = c(0,max(known)*1.3),
             ylab = 'N of genes',
             cex.axis = 0.8)

#text
text(bar,known+max(known)/15,labels = known,cex = 0.8)



#additional
additional<-as.numeric(table(df$group))[c(4,3)]

bar<-barplot(additional, border=F, 
             names.arg=c('PrediXcan','XT-SCAN'),cex.names = 0.8,
             las=2, 
             col=c(light_col['xt'],light_col['st']), 
             main=paste0("Additional \n genes"),cex.main=1,
             ylim = c(0,max(additional)*1.3),
             ylab = 'N of genes',
             cex.axis = 0.8)

#text
text(bar,additional+max(additional)/15,labels = additional,cex = 0.8)

dev.off()


