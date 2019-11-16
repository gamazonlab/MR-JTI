
#single tissue predixcan
#gtex v8

args<-as.numeric(commandArgs(TRUE))
library('glmnet')
folder='st' #predixcan

#set up sub jobs
run_id<-1
run_list<-list()
for (i in 1:49){  #tissues
  for (j in 1:60){  #subjobs per tissue
    run_list[[run_id]]<-c(i,j)
    run_id=run_id+1
  }
}

run_i<-run_list[[args[1]]]
print(run_i)

#tissue
tissue_list<-dir('/data/coxvgi/zhoud2/data/gtex/exp/v8/weighted/')
tissue<-tissue_list[run_i[1]]

#mkdir
com<-paste0('mkdir /data/coxvgi/zhoud2/projects/gtex/weights/',folder,'/ ; mkdir /data/coxvgi/zhoud2/projects/gtex/weights/',folder,'/',tissue)
system(command = com, wait = T)

#get gene list
exp_list<-sub('....$','',dir(paste0('/data/coxvgi/zhoud2/data/gtex/exp/v8/weighted/',tissue)))
geno_list<-sub('....$','',dir('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/dosage_1m/'))
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
    geno<-readRDS(paste0('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/dosage_1m/',gene_list[i],'.rds'))
    geno[,1]<-paste0('GTEX.',geno[,1])
    snp_info<-readRDS(paste0('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/info_1m/',gene_list[i],'.rds'))
    
    #input expression levels
    exp<-readRDS(paste0('/data/coxvgi/zhoud2/data/gtex/exp/v8/weighted/',tissue,'/',gene_list[i],'.rds'))
    exp<-exp[exp$tissue==tissue,] #only keeps expression in target tissue
    
    #merge
    d<-merge(exp,geno,by='sampleid')
    
    #-------------------
    if(ncol(d)>5){ #skip genes with only one SNP
      y<-as.matrix(d[,2]); x<-as.matrix(d[,c(5:ncol(d))]) #assign x and y
      set.seed(as.numeric(sub('^....','',gene_list[i]))) #use ENSG id as seed
      #elastic net
      fit<-cv.glmnet(x=x,y=y, nfolds = 5,keep = T,alpha=0.5,nlambda=50,pmax=200)  
      fit.df <- data.frame(fit$cvm,fit$lambda,1:length(fit$cvm))
      best.lam <- fit.df[which.min(fit.df[,1]),] #find lambda with min cv error
      nrow.best = best.lam[,3] #position of best lambda in cv.glmnet output
      ret <- as.numeric(fit$glmnet.fit$beta[,nrow.best]) #beta
      lambda<-as.numeric(fit$lambda.min) #lambda
    }else{next}
    
    #correlation of predicted and observed expression
    pred_exp = fit$fit.preval[,which(fit$lambda == fit$lambda.min)] #predicted expression
    cor_test_t<-cor.test(y,pred_exp) #correlation test
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
    if (nrow(output)>0 & p<0.05 & r>0.1){ #only output 'imputable Genes'
      out_path<-paste0('/data/coxvgi/zhoud2/projects/gtex/weights/',folder,'/',tissue,'/',gene_list[i],'.rds')
      saveRDS(output,file=out_path)
    }
  }
}else{
  print('out of range')
}







