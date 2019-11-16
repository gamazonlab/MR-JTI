
#single tissue predixcan
#gtex v8

args<-as.numeric(commandArgs(TRUE))
library('glmnet')
main_dir<-'/data/coxvgi/zhoud2/projects/cross_tissue/'
v8_dir<-'/data/coxvgi/zhoud2/projects/v8/'

run_id<-1
run_list<-list()
for (i in 1:49){ #tissue id
  for (j in 1:50){ #subjobs for each tissue
    run_list[[run_id]]<-c(i,j)
    run_id=run_id+1
  }
}

run_i<-run_list[[args[1]]]
print(run_i)

folder='st' #predixcan

#tissue name
tissue_list<-dir(paste0(main_dir,'exp_residual/'))
tissue_list<-as.character(sapply(tissue_list,function(x) strsplit(x,"[.]")[[1]][1]))
tissue<-tissue_list[run_i[1]]

#mkdir for output
com<-paste0('mkdir ',v8_dir,'weights_cv/',folder,'/ ; mkdir ',v8_dir,'weights_cv/',folder,'/',tissue)
system(command = com, wait = T)

#get gene list
exp_list<-sub('....$','',dir(paste0(main_dir,'exp_v8/',tissue)))  #expression gene list
geno_list<-dir('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_1m/')  #genotype gene list
gene_list<-intersect(exp_list,geno_list)

#subjob start and end id
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
    geno<-read.table(paste0('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_1m/',gene_list[i]),header=T,stringsAsFactors =F)
    geno[,1]<-paste0('GTEX.',geno[,1])
    snp_info<-read.table(paste0('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_1m_info/',gene_list[i]),header=T,stringsAsFactors =F)
    
    #input expression levels
    exp<-readRDS(paste0(main_dir,'exp_v8/',tissue,'/',gene_list[i],'.rds'))
    exp<-exp[exp$tissue==tissue,]
    
    #merge
    d<-merge(exp,geno,by='sampleid')
    
    #-------------------
    if(ncol(d)>5){  #skip genes with only one SNP.
      y<-as.matrix(d[,2]); x<-as.matrix(d[,c(5:ncol(d))]) #assign x and y
      set.seed(as.numeric(sub('^....','',gene_list[i]))) #use ensembl id as seed
      #elastic net model training using 5-fold cv
      fit<-cv.glmnet(x=x,y=y, nfolds = 5,keep = T,alpha=0.5,nlambda=50,pmax=200) 
      fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm))
      best.lam <- fit.df[which.min(fit.df[,1]),] #find lambda with min cv error
      nrow.best = best.lam[,3] #position of best lambda in cv.glmnet output
      ret <- as.numeric(fit$glmnet.fit$beta[,nrow.best]) #beta
      lambda<-as.numeric(fit$lambda.min) #lambda
    }else{next}
    
    #correlation of predicted and observed expression
    pred_exp = fit$fit.preval[,which(fit$lambda == fit$lambda.min)] #predicted expression
    cor_test_t<-cor.test(y,pred_exp)  #correlation test
    r<-as.numeric(cor_test_t$estimate)
    p<-as.numeric(cor_test_t$p.value)
    
    #output dataframe
    snp_info$weight<-ret
    snp_info$gene<-gene_list[i]
    snp_info$r2<-r^2
    snp_info$p<-p
    snp_info$lambda=lambda
    snp_info<-snp_info[,c(6,1,2,3,4,5,7,8,9)]
    snp_info<-snp_info[snp_info$weight!=0,]
    output<-snp_info
    
    #writeRDS
    if (nrow(output)>0 & p<0.05 & r>0.1){ #only output iGene defined as p<0.05 $ r>0.10
      out_path<-paste0(v8_dir,'weights_cv/',folder,'/',tissue,'/',gene_list[i],'.rds')
      saveRDS(output,file=out_path)
    }
  }
}else{
  print('out of range')
}







