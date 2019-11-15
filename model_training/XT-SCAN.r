
#cross tissue XT-SCAN
#gtex v8

args<-as.numeric(commandArgs(TRUE))
library('glmnet')

main_dir<-'/data/coxvgi/zhoud2/projects/cross_tissue/'
v8_dir<-'/data/coxvgi/zhoud2/projects/v8/'

run_id<-1
run_list<-list()
for (i in 1:49){ #1:49
  for (j in 1:400){  #50
    run_list[[run_id]]<-c(i,j)
    run_id=run_id+1
  }
}

run_i<-run_list[[args[1]]]
print(run_i)

folder='xt_cv'

#tissue
tissue_list<-dir(paste0(main_dir,'exp_residual/'))
tissue_list<-as.character(sapply(tissue_list,function(x) strsplit(x,"[.]")[[1]][1]))
tissue<-tissue_list[run_i[1]]

#mkdir
dir.create(paste0(v8_dir,'weights_cv/',folder,'/'))
dir.create(paste0(v8_dir,'weights_cv/',folder,'/',tissue))

#load tissue level dhs similarity
dhs_weight_tissue<-read.table(paste0('/data/coxvgi/zhoud2/projects/cross_tissue/encode/simularity_tissue/',tissue),header = T,stringsAsFactors = F)
#get gene level dhs similarity gene list
dhs_gene_list<-dir(paste0('/data/coxvgi/zhoud2/projects/cross_tissue/encode/simularity/',tissue,'/gene/'))

#get gene list
exp_list<-sub('....$','',dir(paste0(main_dir,'exp_v8/',tissue))) #expression 
geno_list<-dir('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_1m/') #dosage genotype
gene_list<-intersect(exp_list,geno_list)

#start_id
i_start=(run_i[2]-1)*50+1  
i_end=run_i[2]*50  
if (i_end>length(gene_list)){
  i_end=length(gene_list)
}

#---------------------------
#i=which(gene_list=='ENSG00000240654')

if (i_start<length(gene_list)){
  for (i in i_start:i_end){
    print(paste0('INFO processing ',gene_list[i],' ...'))
    
    #---data prepare---
    
    #load geno snpinfo
    #1m
    geno_1m<-read.table(paste0('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_1m/',gene_list[i]),header=T,stringsAsFactors =F)
    geno_1m[,1]<-paste0('GTEX.',geno_1m[,1])
    snp_info_1m<-read.table(paste0('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_1m_info/',gene_list[i]),header=T,stringsAsFactors =F)
    if(ncol(geno_1m)==2){next}
    #100k
    geno_100k<-try(read.table(paste0('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_100k/',gene_list[i]),header=T,stringsAsFactors =F))
    if('try-error' %in% class(geno_100k)){
      geno_100k=geno_1m
      snp_info_100k=snp_info_1m
    }else if(ncol(geno_100k)==2){
      geno_100k=geno_1m
      snp_info_100k=snp_info_1m
    }else{
      geno_100k[,1]<-paste0('GTEX.',geno_100k[,1])
      snp_info_100k<-read.table(paste0('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_100k_info/',gene_list[i]),header=T,stringsAsFactors =F)
    }
    
    #load expression levels
    exp<-readRDS(paste0(main_dir,'exp_v8/',tissue,'/',gene_list[i],'.rds'));colnames(exp)[3]<-'exp_w' 
    
    #merge with dhs weights
    if (gene_list[i] %in% dhs_gene_list){
      dhs_weight_gene<-read.table(paste0('/data/coxvgi/zhoud2/projects/cross_tissue/encode/simularity/',tissue,'/gene/',gene_list[i]),header = T,stringsAsFactors = F)
      if (is.na(dhs_weight_gene[1,2])){
        dhs_weight<-dhs_weight_tissue
      }else{
        dhs_weight_gene[which(dhs_weight_gene[,2]<0),2]<-0
        colnames(dhs_weight_gene)<-c('tissue','dhs_w')
        dhs_weight<-dhs_weight_gene
      }
    }else{
      dhs_weight<-dhs_weight_tissue
    }
    exp<-merge(exp,dhs_weight,by='tissue')  
    #col: tissue sampleid exp exp_w dhs_w
    
    
    #---model training---
    
    #fit single tissue model to get proper window size and lambda range for each gene
    exp_st<-exp[exp$tissue==tissue,] #load single tissue data
    d_st_100k<-merge(exp_st,geno_100k,by='sampleid')
    set.seed(as.numeric(sub('^....','',gene_list[i])))
    fit_100k<-cv.glmnet(x=as.matrix(d_st_100k[,6:ncol(d_st_100k)]),y=as.matrix(d_st_100k[,'exp']), nfolds = 5,keep = T,alpha=0.5,nlambda=50,pmax=200)
    
    d_st_1m<-merge(exp_st,geno_1m,by='sampleid')
    set.seed(as.numeric(sub('^....','',gene_list[i])))
    fit_1m<-cv.glmnet(x=as.matrix(d_st_1m[,6:ncol(d_st_1m)]),y=as.matrix(d_st_1m[,'exp']), nfolds = 5,keep = T,alpha=0.5,nlambda=50,pmax=200)
    
    #find proper window size and set up for downstream analysis
    #The cross tissue lambda range was set around the best single tissue lambda
    if(which.min(c(min(fit_100k$cvm),min(fit_1m$cvm)))==1){
      windows='100k';geno<-geno_100k;snp_info<-snp_info_100k
      lambda_list<-fit_100k$lambda[max(1,which(fit_100k$lambda==fit_100k$lambda.min)-10):min(length(fit_100k$lambda),which(fit_100k$lambda==fit_100k$lambda.min)+10)]
    }else{
      windows='1m';geno<-geno_1m;snp_info<-snp_info_1m
      lambda_list<-fit_1m$lambda[max(1,which(fit_1m$lambda==fit_1m$lambda.min)-10):min(length(fit_1m$lambda),which(fit_1m$lambda==fit_1m$lambda.min)+10)]
    }
    
    #dataframe for grid search and performance collection
    r_matrix<-as.data.frame(matrix(data=NA,nrow=(4*4),ncol=6))
    i_loop<-1
    for (exp_power in c(1,4,16,64)){
      for (dhs_power in c(1,4,16,64)){
        r_matrix[i_loop,1:2]<-c(exp_power,dhs_power)
        i_loop=i_loop+1
      }
    }
    colnames(r_matrix)<-c('exp_power','dhs_power','lambda','window','r_test','p_test')
    r_matrix$window=windows
    
    
    #sample id map for each fold (fold=5)
    d_tmp<-merge(exp,geno_100k,by='sampleid')
    sample_all<-unique(d_tmp$sampleid)
    set.seed(as.numeric(sub('^....','',gene_list[1]))+1)
    id_map<-data.frame(sampleid=sample_all,id_map=sample(rep(seq(1,5),ceiling(length(sample_all)/5)),length(sample_all)),stringsAsFactors = F)
    
    #beta list of each hyper-parameter pairs
    beta_list<-list()
    
    #go through each of the hyper-parameter pairs
    for (j in 1:nrow(r_matrix)){
      print(j)
      exp_power=r_matrix[j,1];dhs_power=r_matrix[j,2]
      exp$w<-(exp$exp_w)^exp_power*(exp$dhs_w)^dhs_power
      
      #merge
      d<-merge(exp,geno,by='sampleid')
      
      #rm tissues with w<0.1
      d<-d[d[,'w']>=0.1,]; d<-d[!is.na(d[,1]),]
      
      #id map
      d<-merge(id_map,d,by='sampleid')
      
      #target tissue position (the performance will be only estimated in the target tissue)
      tt_pos<-which(d$tissue==tissue)
      
      #cross-tissue weighted elastic net
      set.seed(as.numeric(sub('^....','',gene_list[1]))+2)
      ans<-cv.glmnet(x=as.matrix(d[,8:ncol(d)]),y=as.matrix(d[,'exp']),weights=d[,'w'],foldid=d[,'id_map'],pmax=200,lambda=lambda_list,keep = T)
      
      #correlation between pred and obs expression levels
      cor_ans<-cor.test(ans$fit.preval[tt_pos,which(ans$lambda==ans$lambda.min)],d[tt_pos,'exp'])
      
      r_matrix[j,'r_test']<-cor_ans$estimate
      r_matrix[j,'p_test']<-cor_ans$p.value
      r_matrix[j,'lambda']<-ans$lambda.min
      
      #collect the weights
      beta_list[[j]]<-as.numeric(ans$glmnet.fit$beta[,which(ans$lambda==ans$lambda.min)])
    }
    
    #find the best hyperparameters with the largest r_test
    best_row=which.max(r_matrix[,'r_test'])[1]
    r_test=r_matrix[best_row,'r_test']
    p_test=r_matrix[best_row,'p_test']
    
    #---output---
    
    if (p_test<0.05 & r_test>0.1){
      snp_info$gene<-gene_list[i]
      snp_info$r2<-r_test^2
      snp_info$p<-p_test
      snp_info$lambda<-r_matrix[best_row,'lambda']
      snp_info$weight<-beta_list[[best_row]]
      snp_info<-snp_info[,c(5,1:4,9,6:8)]
      nocv<-snp_info[snp_info$weight!=0,]
      if (nrow(nocv)>0){
        out_path<-paste0(v8_dir,'weights_cv/',folder,'/',tissue,'/',gene_list[i],'.rds') 
        saveRDS(nocv,file=out_path)
      }
    }
    
  }
}else{
  print('out of range')
}









