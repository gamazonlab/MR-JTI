args=as.numeric(commandArgs(TRUE)) 
library(ggplot2)

tissue_list<-dir('/data/c***/z***/projects/gtex/weights/st/')
tissue=tissue_list[args]


main_path<-'/data/c***/z***/data/biovu/pred_exp_23k/'
out_path<-'/data/c***/z***/projects/gtex/simulation/type1error_genome/'
simu_times=1  #test 50 100 1000

#genelist
gene_list_st<-dir(paste0(main_path,'st/',tissue)) #PrediXcan iGenes
gene_list_ct<-dir(paste0(main_path,'xt/',tissue)) #JTI iGenes
gene_list_ut<-dir(paste0(main_path,'ut/',tissue)) #UTMOST iGenes

#creat dataframes to collect results
out_s<-out_c<-out_u<-as.data.frame(matrix(data=NA,ncol=3,nrow=0))
colnames(out_s)<-colnames(out_c)<-colnames(out_u)<-c('gene','i','p')

#st
for (i in 1:length(gene_list_st)){
  #for (i in 1:100){
  print(i)
  
  gene<-gene_list_st[i]
  #load predicted expression
  d_s<-readRDS(paste0(main_path,'st/',tissue,'/',gene)) #single tissue
  
  for (j in 1:simu_times){   #1000 times
    #print(paste0(i,' ',j))
    
    #single tissue
    y=rnorm(nrow(d_s),mean=0,sd=1) #e~N(0,1^2)
    fit_s<-summary(lm(y~d_s[,2]))
    out_s[j,2]<-j 
    out_s[j,3]<-fit_s$coefficients[2,4]  #p value
    
  }
  
  out_s[,1]<-gene
  
  #combine results
  if (i==1){
    out_s_o<-out_s
    
  }else{
    out_s_o<-rbind(out_s_o,out_s)
  }
}




#ct
for (i in 1:length(gene_list_ct)){
  #for (i in 1:100){
  print(i)
  
  gene<-gene_list_ct[i]
  #load predicted expression
  d_c<-readRDS(paste0(main_path,'xt/',tissue,'/',gene)) #single tissue
  
  for (j in 1:simu_times){   #1000 times
    #print(paste0(i,' ',j))
    if(var(d_c[,2])==0){next}
    #single tissue
    y=rnorm(nrow(d_c),mean=0,sd=1) #e~N(0,1^2)
    fit_c<-summary(lm(y~d_c[,2]))
    out_c[j,2]<-j 
    out_c[j,3]<-fit_c$coefficients[2,4]  #p value
    
  }
  
  out_c[,1]<-gene
  
  #combine results
  if (i==1){
    out_c_o<-out_c
    
  }else{
    out_c_o<-rbind(out_c_o,out_c)
  }
}

#ut
for (i in 1:length(gene_list_ut)){
  #for (i in 1:100){
  print(i)
  
  gene<-gene_list_ut[i]
  #load predicted expression
  d_u<-readRDS(paste0(main_path,'ut/',tissue,'/',gene)) #single tissue
  
  for (j in 1:simu_times){   #1000 times
    #print(paste0(i,' ',j))
    if(var(d_u[,2])==0){next}
    #single tissue
    y=rnorm(nrow(d_u),mean=0,sd=1) #e~N(0,1^2)
    fit_u<-summary(lm(y~d_u[,2]))
    out_u[j,2]<-j 
    out_u[j,3]<-fit_u$coefficients[2,4]  #p value
    
  }
  
  out_u[,1]<-gene
  
  #combine results
  if (i==1){
    out_u_o<-out_u
    
  }else{
    out_u_o<-rbind(out_u_o,out_u)
  }
}




out_s_o<-out_s_o[!is.na(out_s_o$p),]

out_c_o<-out_c_o[!is.na(out_c_o$p),]

out_u_o<-out_u_o[!is.na(out_u_o$p),]




gg_qqplot <- function(ps, ci = 0.95,title='QQ-plot',ymax) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 1.5) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2,color='dodgerblue3') +
    geom_line(aes(expected, clower), linetype = 2,color='dodgerblue3') +
    ylim(0,ymax) +
    xlab(log10Pe) +
    ylab(log10Po) +
    ggtitle(title)
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


for (i in 1:10){
  tissue<-sub('_',' ',tissue)
}


png(paste0(out_path,'plot/',tissue,'_allgenes.png'),width = 2800,height = 1000,res=200)
#par(mfcol=c(1,1))

#out_s_o$p<-10^(out_s_o$p*-1)
#out_c_o$p<-10^(out_c_o$p*-1)


p1=gg_qqplot(out_s_o$p,ymax=6,title=paste0(tissue,' (PrediXcan)'))
p2=gg_qqplot(out_u_o$p,ymax=6,title=paste0(tissue,' (UTMOST)'))
p3=gg_qqplot(out_c_o$p,ymax=6,title=paste0(tissue,' (JTI)'))


multiplot(p1, p2, p3, cols=3)

dev.off()

#--output results--

write.table(out_s_o,paste0(out_path,'st_',tissue,'_allgenes.txt'),sep='\t',quote = F,row.names = F)
write.table(out_c_o,paste0(out_path,'xt_',tissue,'_allgenes.txt'),sep='\t',quote = F,row.names = F)
write.table(out_u_o,paste0(out_path,'ut_',tissue,'_allgenes.txt'),sep='\t',quote = F,row.names = F)







