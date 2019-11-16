
require("RColorBrewer")
require('directlabels')
library(reshape2)
library(knitr)
options(stringsAsFactors = F)
#plot mr result

df<-read.csv("~/../Dropbox/DansPaper/data/mr/mr.csv",header = T,stringsAsFactors = F)
df<-df[,c('geneid','sig_real_ci')]

#---load annotation---
gene_anno_pos<-read.table('~/../Dropbox/annotation/encode/gencode.v27.gene',header = F);gene_anno_pos<-gene_anno_pos[gene_anno_pos$V1 %in% c(paste0('chr',seq(1,22))),];gene_anno_pos$pos<-(gene_anno_pos$V3+gene_anno_pos$V4)/2;gene_anno_pos$chr<-as.numeric(sub('^...','',gene_anno_pos[,1]));gene_anno_pos$geneid<-sapply(gene_anno_pos$V5,function(x) strsplit(x,"[.]")[[1]][1]);gene_anno_pos<-gene_anno_pos[,c('geneid','chr','pos')]
gene_anno_name<-read.table('~/../Dropbox/annotation/encode/v27.anno.gene',header = F);gene_anno_name$geneid<-sapply(gene_anno_name$V2,function(x) strsplit(x,"[.]")[[1]][1]);gene_anno_name$genename<-gene_anno_name$V6;gene_anno_name<-gene_anno_name[,c('geneid','genename')]

gene_anno<-merge(gene_anno_name,gene_anno_pos,by='geneid')

#---merge annotation---
df<-merge(gene_anno,df,by=1,all=T)
df$pos_plot=df$pos

#---generate_pos---
chr_box<-as.data.frame(matrix(data=NA,nrow=22,ncol=3))
chr_df<-gene_anno[gene_anno$chr==1,]
chr_box[1,1]<-max(chr_df$pos)
chr_box[1,2]<-0
chr_box[1,3]<-chr_box[1,1]/2

for (i in 2:22){
  chr_df<-gene_anno[gene_anno$chr==i,]
  chr_box[i,1]<-max(chr_df$pos)
  chr_box[i,2]<-max(chr_df$pos)+chr_box[i-1,2]
  df[which(df$chr==i),'pos_plot']<-df[which(df$chr==i),'pos_plot']+chr_box[i,2]
  chr_box[i,3]<-chr_box[i,2]+chr_box[i,1]/2
}



#---rename cols---
colnames(df)<-c('geneid','genename','chr','pos_old','sig','pos')

#--load twas--
twas<-read.csv('~/../Dropbox/DansPaper/data/apply_gwas/ukbb_xt_Liver.csv')

pvalue.extreme <- function(z) {
  log.pvalue <- log(2) + pnorm(abs(z), lower.tail = FALSE, log.p = TRUE)
  log10.pvalue <- log.pvalue/log(10) ## from natural log to log10
  mantissa <- 10^(log10.pvalue %% 1)
  exponent <- log10.pvalue %/% 1
  
  logp<-log10.pvalue*-1
  return(logp)
  ##return(sprintf("p value is %1.2f times 10^(%d)",mantissa,exponent))
}
twas$logp<-pvalue.extreme(twas$zscore)

twas<-twas[,c(1,4,5,13)]
colnames(twas)<-c('geneid','beta_twas','p_twas','logp')

#--merge with twas--
df<-merge(df,twas,by='geneid')
df<-df[,c('geneid','genename','chr','pos_old','pos','beta_twas','p_twas','logp','sig')]

#---load know---
know<-read.csv('~/../Dropbox/XTSCAN/glgc/knowgenes.csv',header = T,stringsAsFactors = F)


#--gene-group--
df$group<-NA
df[df$p_twas>=0.05,'group']<-'TWAS non sig'
df[which(df$p_twas<0.05 & df$sig!='sig'),'group']<-'TWAS sig (only)'
df[which(df$p_twas<0.05 & is.na(df$sig)),'group']<-'TWAS sig (only)' #'gold3'
df[which(df$p_twas<0.05 & df$sig=='sig'),'group']<-'MR & TWAS sig'
df$genename_highlight<-ifelse((df$group=='MR & TWAS sig' & df$genename %in% c(know$genename)),df$genename,'ZZZ')
#'UQCC','TAGLN','CSNK1G3','LSM7','POLK','CELSR2','PSRC1','TDRD15'

df<-df[df$p_twas<0.50,]

#reorder legend
df$group<-as.factor(df$group)
df$group <- factor(df$group, levels = c("MR & TWAS sig", "TWAS sig (only)", "TWAS non sig"))

#plot 
library(ggrepel)

#order
df<-df[order(df$genename_highlight,decreasing = T),]
df[,'genename_highlight']<-ifelse(df[,'genename_highlight']=='ZZZ',NA,df[,'genename_highlight'])


png('~/../Dropbox/DansPaper/figures/tmp/mr_manha.png',width = 1800,height = 1000,res=250)

set.seed(1)
ggplot(df, aes(x=pos, y=logp, label=genename_highlight)) +
  ylim(0,200) +
  scale_x_continuous(breaks=chr_box$V3, labels = c(as.character(seq(1,15)),' ','17','','19','','21','')) +  #x axis
  scale_y_continuous(breaks=c(1,2,3,5,10,20,30,50,100,200,300,500),trans='log10') +  #y axis
  geom_point(data = df, aes(x = pos, y = logp, color = group, fill=group, shape=group)) +
  scale_fill_manual(values = c("TWAS non sig"='grey', "MR & TWAS sig"='gold', 'TWAS sig (only)'=brewer.pal(12, "Paired")[4])) +
  scale_color_manual(values = c("TWAS non sig"='grey', "MR & TWAS sig"='black', 'TWAS sig (only)'='black')) +
  scale_shape_manual(values = c("TWAS non sig"=16, "MR & TWAS sig"=21, 'TWAS sig (only)'=21)) +
  
  labs(x = "Chromosome", y = "-Log(p)", title = "") +
  theme_bw() +  #rm background
  theme(panel.grid =element_blank()) +  #rm grids
  theme(legend.position= "right",legend.title = element_blank()) + #set legend
  
  geom_label_repel( #non overlapped labels
    size=2, 
    #col='white',
    fill=rgb(255, 255, 255, 200, maxColorValue=255),
    #nudge_y=200 - df$logp[!is.na(df$genename_highlight)],
    #direction='x', #only work on Y-axis
    #nudge_x=1.5e8, #shift to the right
    segment.alpha = 0.5,  #transparent of segment
    #min.segment.length = 1e9
    #parse=T,
    box.padding=0.5,
    nudge_x=2e8,
    fontface='italic'
  ) 

dev.off()







#-----------bar plot-------------

png('~/../Dropbox/DansPaper/figures/tmp/mr_bar.png',width = 2000,height = 1000,res=400)
par(mfcol=c(1,4),mar=c(3,4,3,1),mgp=c(2.5,0.75,0.25))

#overall
overall<-c(length(which(df$group=='TWAS non sig')),length(which(df$group=='TWAS sig (only)')),length(which(df$group=='MR & TWAS sig')))

bar<-barplot(overall, border=F, 
             #names.arg=c('TWAS non sig','TWAS sig (only)','TWAS & MR sig'),
             las=2, 
             col=c('grey', brewer.pal(12, "Paired")[4],'gold'), 
             main=paste0("Overall"),cex.main=0.8,
             ylim = c(0,max(overall)*1.2),
             ylab = 'N of genes',
             cex.axis = 0.8)
#label x
text(bar,seq(1,3)+0.5, srt = 30, adj= 1, xpd = TRUE, labels = c('TWAS non sig ','TWAS sig (only) ','TWAS & MR sig '), cex=0.8)
#text
text(bar,overall+100,labels = overall,cex = 0.8)



#known / additional
know_gene<-know[,'geneid']
twas_non_sig_gene<-df[df$group=='TWAS non sig','geneid']
twas_sig_gene<-df[df$group=='TWAS sig (only)','geneid']
mr_sig_gene<-df[df$group=='MR & TWAS sig','geneid']

#twas non sig
twas_non_sig<-c(length(intersect(know_gene,twas_non_sig_gene)),length(twas_non_sig_gene)-length(intersect(know_gene,twas_non_sig_gene)))

bar<-barplot(twas_non_sig , border=F , 
             #names.arg=c('TWAS non sig','TWAS sig (only)','TWAS & MR sig'),
             las=2, 
             col='grey', 
             main=paste0("TWAS non sig \n ",round(twas_non_sig[1]/(twas_non_sig[1]+twas_non_sig[2])*100,2),'% known genes'),cex.main=0.8,
             ylim = c(0,max(twas_non_sig)*1.2),
             ylab = 'N of genes',
             cex.axis = 0.8)
#label x
text(bar,seq(1,2)+0.5, srt = 30, adj= 1, xpd = TRUE, labels = c('Known ','Others '), cex=0.8)
#text
text(bar,twas_non_sig+max(twas_non_sig)/15,labels = twas_non_sig,cex = 0.8)

#twas sig
twas_sig<-c(length(intersect(know_gene,twas_sig_gene)),length(twas_sig_gene)-length(intersect(know_gene,twas_sig_gene)))

bar<-barplot(twas_sig , border=F , 
             #names.arg=c('TWAS non sig','TWAS sig (only)','TWAS & MR sig'),
             las=2, 
             col=brewer.pal(12, "Paired")[4], 
             main=paste0("TWAS sig (only) \n ",round(twas_sig[1]/(twas_sig[1]+twas_sig[2])*100,2),'% known genes'),cex.main=0.8,
             ylim = c(0,max(twas_sig)*1.2),
             ylab = 'N of genes',
             cex.axis = 0.8)
#label x
text(bar,seq(1,2)+0.5, srt = 30, adj= 1, xpd = TRUE, labels = c('Known ','Others '), cex=0.8)
#text
text(bar,twas_sig+max(twas_sig)/15,labels = twas_sig,cex = 0.8)


#mr sig
mr_sig<-c(length(intersect(know_gene,mr_sig_gene)),length(mr_sig_gene)-length(intersect(know_gene,mr_sig_gene)))

bar<-barplot(mr_sig , border=F , 
             #names.arg=c('TWAS non sig','TWAS sig (only)','TWAS & MR sig'),
             las=2, 
             col='gold', 
             main=paste0("TWAS & MR sig \n ",round(mr_sig[1]/(mr_sig[1]+mr_sig[2])*100,2),'% known genes'),cex.main=0.8,
             ylim = c(0,max(mr_sig)*1.2),
             ylab = 'N of genes',
             cex.axis = 0.8)
#label x
text(bar,seq(1,2)+0.5, srt = 30, adj= 1, xpd = TRUE, labels = c('Known ','Others '), cex=0.8)
#text
text(bar,mr_sig+max(mr_sig)/15,labels = mr_sig,cex = 0.8)



dev.off()































