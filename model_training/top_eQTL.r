
#single tissue top
#gtex v8

args<-as.numeric(commandArgs(TRUE))

folder='top'

#set up sub jobs
run_id<-1
run_list<-list()
for (i in c(22)){  #tissues
  for (j in 1:40){  #subjobs per tissue
    run_list[[run_id]]<-c(i,j)
    run_id=run_id+1
  }
}

run_i<-run_list[[args[1]]]
print(run_i)

#tissue
tissue_list<-dir('/data/c***/z***/data/gtex/exp/v8/weighted/')
tissue<-tissue_list[run_i[1]]

#mkdir
com<-paste0('mkdir /data/c***/z***/projects/gtex/weights/',folder,'/ ; mkdir /data/c***/z***/projects/gtex/weights/',folder,'/',tissue)
system(command = com, wait = T)

#get gene list
exp_list<-sub('....$','',dir(paste0('/data/c***/z***/data/gtex/exp/v8/weighted/',tissue)))
geno_list<-sub('....$','',dir('/data/c***/z***/data/gtex/geno/v8/gene/dosage_1m/'))
gene_list<-intersect(exp_list,geno_list)

#start and end id
i_start=(run_i[2]-1)*500+1
i_end=run_i[2]*500
if (i_end>length(gene_list)){
  i_end=length(gene_list)
}

#---------------------------

output<-as.data.frame(matrix(data=NA,nrow = 0,ncol=8))

if (i_start<length(gene_list)){
  for (i in i_start:i_end){
    print(i)
    #input geno and snp info
    geno<-readRDS(paste0('/data/c***/z***/data/gtex/geno/v8/gene/dosage_1m/',gene_list[i],'.rds'))
    geno[,1]<-paste0('GTEX.',geno[,1])
    snp_info<-readRDS(paste0('/data/c***/z***/data/gtex/geno/v8/gene/info_1m/',gene_list[i],'.rds'))
    
    #input expression levels
    exp<-readRDS(paste0('/data/c***/z***/data/gtex/exp/v8/weighted/',tissue,'/',gene_list[i],'.rds'))
    exp<-exp[exp$tissue==tissue,] #only keeps expression in target tissue
    
    #merge
    d<-merge(exp,geno,by='sampleid')
    
    
    #find the top
    find_top_p<-c()
    find_top_beta<-c()
    find_top_r<-c()
    
    for (j in 1:(ncol(d)-4)){
      ans<-summary(lm(d$exp~d[,j+4]))
      find_top_beta[j]<-ans$coefficients['d[, j + 4]','Estimate']
      find_top_p[j]<-ans$coefficients['d[, j + 4]','Pr(>|t|)']
      find_top_r[j]<-abs(cor(d$exp,d[,j+4]))
    }
    
    top<-which.min(find_top_p)
    rsid<-colnames(d)[top+4]
    
    snp_info<-snp_info[snp_info$rsid==rsid,]
    snp_info$weight<-find_top_beta[top]
    snp_info$gene<-gene_list[i]
    snp_info$r2<-find_top_r[top]^2
    snp_info$p<-find_top_p[top]
    snp_info$lambda=0
    snp_info<-snp_info[,c(6,1,2,3,4,5,7,8,9)]
    snp_info<-snp_info[snp_info$weight!=0,]
    output<-snp_info
    
    #writeRDS
    if (nrow(output)>0 & find_top_p[top]<0.05 & find_top_r[top]>0.1){
      out_path<-paste0('/data/c***/z***/projects/gtex/weights/',folder,'/',tissue,'/',gene_list[i],'.rds')
      saveRDS(output,file=out_path)
    }
  }
}else{
  print('out of range')
}








