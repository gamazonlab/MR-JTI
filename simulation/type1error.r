args=as.numeric(commandArgs(TRUE)) 
library(ggplot2)

tissue_list<-dir('/data/coxvgi/zhoud2/projects/gtex/weights/st/')
tissue=tissue_list[args]

main_path<-'/data/coxvgi/zhoud2/data/biovu/pred_exp_23k/' #BioVU real genotyping data
out_path<-'/data/coxvgi/zhoud2/projects/gtex/simulation/type1error/'
simu_times=1000  #test 50 100 1000

#genelist
gene_list_st<-dir(paste0(main_path,'st/',tissue)) #PrediXcan iGenes
gene_list_ct<-dir(paste0(main_path,'xt/',tissue)) #XT-SCAN iGenes
gene_list_all<-intersect(gene_list_ct,gene_list_st)

#100 random genes
set.seed(args)
gene_list<-sample(gene_list_all,100,replace = F)

#creat dataframes to collect results
out_s<-out_c<-as.data.frame(matrix(data=NA,ncol=4,nrow=0))
colnames(out_s)<-colnames(out_c)<-c('gene','i','p','exp_p')

for (i in 1:length(gene_list)){
  print(i)
  
  gene<-gene_list[i]
  #load predicted expression
  d_s<-readRDS(paste0(main_path,'st/',tissue,'/',gene)) #single tissue
  d_c<-readRDS(paste0(main_path,'xt/',tissue,'/',gene)) #corss tissue
  
  #expected -logp
  exp_p<-c()
  for (k in 1:simu_times){
    y_r=rnorm(nrow(d_s),mean=0,sd=1)
    x_r=rnorm(nrow(d_s),mean=0,sd=1)
    fit_exp<-summary(lm(y_r~x_r))
    exp_p[k]<-fit_exp$coefficients[2,4]
  }
  exp_p<-sort(exp_p)
  
  
  for (j in 1:simu_times){   #1000 times
    #print(paste0(i,' ',j))
    
    #single tissue
    y=rnorm(nrow(d_s),mean=0,sd=1) #e~N(0,1^2)
    fit_s<-summary(lm(y~d_s[,2]))
    out_s[j,2]<-j 
    out_s[j,3]<-fit_s$coefficients[2,4]  #p value
    
    #cross tissue
    y=rnorm(nrow(d_c),mean=0,sd=1) #e~N(0,1^2)
    fit_c<-summary(lm(y~d_c[,2]))
    out_c[j,2]<-j
    out_c[j,3]<-fit_c$coefficients[2,4]  #p value
  }
  
  out_s<-out_s[order(out_s[,3]),]
  out_c<-out_c[order(out_c[,3]),]
  
  out_s[,1]<-out_c[,1]<-gene
  out_s[,4]<-out_c[,4]<-exp_p
  
  #--ci--
  ci=0.95
  n=simu_times
  
  st_tmp <- data.frame(
    observed = -log10(sort(out_s$p)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  
  ct_tmp <- data.frame(
    observed = -log10(sort(out_c$p)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  
  #combine results
  if (i==1){
    st<-st_tmp
    ct<-ct_tmp
    out_s_o<-out_s
    out_c_o<-out_c
  }else{
    st<-rbind(st,st_tmp)
    ct<-rbind(ct,ct_tmp)
    out_s_o<-rbind(out_s_o,out_s)
    out_c_o<-rbind(out_c_o,out_c)
  }
  
}

#--qq plot--

gg_qqplot <- function(df, ci = 0.95,title='ggplot') {
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df,dpi=300) +
    geom_point(aes(expected, observed), shape = 1, size = 1) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2,color='dodgerblue3') +
    geom_line(aes(expected, clower), linetype = 2,color='dodgerblue3') +
    xlab(log10Pe) +
    ylab(log10Po) +
    ggtitle(title) +
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"),
          axis.title.x=element_text(size=20,colour = 'black'),
          axis.title.y=element_text(size=20,colour = 'black'),
          axis.text.x=element_text(size=20,colour = 'black'),
          axis.text.y=element_text(size=20,colour = 'black'),
          title=element_text(size=15,colour = 'black'))
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#rm '_' in tissue names
for (i in 1:10){
  tissue<-sub('_',' ',tissue)
}

png(paste0(out_path,'plot/',tissue,'.png'),width = 1000,height = 500)
par(mfcol=c(1,1))

#out_s_o$p<-10^(out_s_o$p*-1)
#out_c_o$p<-10^(out_c_o$p*-1)
p1<-gg_qqplot(st,title=paste0(tissue,' (PrediXcan)'))
p2<-gg_qqplot(ct,title=paste0(tissue,' (XT-SCAN)'))
multiplot(p1, p2, cols=2)

dev.off()

#--output results--

write.table(out_s_o,paste0(out_path,'st_',tissue,'.txt'),sep='\t',quote = F,row.names = F)
write.table(out_c_o,paste0(out_path,'xt_',tissue,'.txt'),sep='\t',quote = F,row.names = F)



















