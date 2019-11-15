
args<-as.numeric(commandArgs(TRUE))  #subjobs

model='xt' #XT-SCAN model
options(digits = 4)
tissue='Liver'
nf=5 #number of fold for cross-validation
prune=0.2  #LD pruning r2 0.2

library("glmnet")
#library("mvtnorm")
library('HDCI')
#library('selectiveInference')
#library('hdi')


#---function Residual Bootstrap LASSO---
RB_LASSO<-function (x, y, B = 500, alpha = 0.05, nfolds = 5, foldid, cv.OLS = FALSE, tau = 0, parallel = FALSE, standardize = TRUE, intercept = TRUE, parallel.boot = FALSE, ncores.boot = 1, ...) {
  # 
  #   B = 500  #no of bootstrap
  #   alpha = 0.05  #95%CI
  #   nfolds = 5  #cv nfolds
  #   cv.OLS = FALSE
  #   tau = 0
  #   parallel = FALSE
  #   standardize = TRUE
  #   intercept = TRUE
  #   parallel.boot = FALSE
  #   ncores.boot = 1
  #   
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  selectset <- rep(0, p)
  Beta <- rep(0, p)
  
  #fit a global model
  globalfit <- glmnet(x, y, standardize = standardize, intercept = intercept)
  
  #get the lambda list
  lambda <- globalfit$lambda
  
  #cv glmnet
  cvfit <- escv.glmnet(x, y, lambda = lambda, nfolds = nfolds, tau = tau, cv.OLS = cv.OLS, parallel = parallel, standardize = standardize,intercept = intercept)
  
  #find the lambda with min cv.error
  lambda.opt <- cvfit$lambda.cv
  
  #beta including intercept
  fitlasso <- predict(globalfit, type = "coefficients",s = lambda.opt)
  #rm the intercept
  Beta <- fitlasso[-1] 
  
  #threshold beta
  Beta<-ifelse(abs(Beta)<1/ncol(x),0,Beta) 
  
  #predicted values
  fit_value <- predict(globalfit, newx = x, s = lambda.opt)
  
  #get residual
  residual <- y - fit_value
  residual_center <- residual - mean(residual)
  
  #beta matrix
  Beta.boot <- matrix(0, nrow = B, ncol = p)
  
  for (i in 1:B) {
    
    #resample bootstrap
    resam <- sample(1:n, n, replace = TRUE)
    ystar <- fit_value + residual_center[resam]
    
    #lasso
    boot.obj <- Lasso(x = x, y = ystar, lambda = lambda.opt, standardize = standardize, intercept = intercept)
    
    #get beta(alpha)
    Beta.boot[i,] <- boot.obj$beta
    
    #threshold beta(alpha)
    Beta.boot[i,]<-ifelse(abs(Beta.boot[i,])<1/ncol(x),0,Beta.boot[i,])
  }
  
  #mean of beta(alpha)
  beta_mean<-apply(Beta.boot,2,function(x) mean(x))
  
  #se of beta(alpha)
  beta_se=apply(Beta.boot, 2, function(x) sqrt(var(x)))
  
  #p value
  pvalue<-pnorm(-abs(beta_mean[1]/beta_se[1]))*2
  
  #95% CI
  interval <- matrix(0, 2, p)
  bound.percentile <- apply(Beta.boot, 2, function(u) {
    quantile(u, prob = c(1 - alpha/2, alpha/2))
  })
  interval[1, ] <- 2 * Beta - bound.percentile[1, ]
  interval[2, ] <- 2 * Beta - bound.percentile[2, ]
  #95% CI significance
  sig <- apply(interval, 2, function(x) ifelse(x[1]*x[2]>0,'sig','nonsig'))
  
  #output
  object <- list(lambda.opt = lambda.opt, beta = beta_mean, interval = interval,pvalue=pvalue)
  object
}


#load LD-pruned snp list 
snp_prune<-read.table(paste0('/gpfs23/data/coxvgi/zhoud2/data/gtex/geno/v8/liftover/pruning/prune',prune,'.prune.in'),stringsAsFactors = F)

#load LDL gwas result
gwas<-read.table('/data/coxvgi/zhoud2/data/ukbb/bn/gwas201807/processed/LDLq.txt',header = T,stringsAsFactors = F)
gwas<-gwas[,c(2,7,3,6)]
gwas[,2]<-toupper(gwas[,2]) #allele
colnames(gwas)<-c('rsid','effect_allele_gwas','beta_gwas','p_gwas')
gwas<-gwas[which(gwas$rsid %in% snp_prune[,1]),]  #all or pruned

#get imputable gene list
gene_list<-sub('....$','',dir(paste0('/data/coxvgi/zhoud2/projects/gtex/weights/',model,'/',tissue)))

#load TWAS results (only run for genes with significant TWAS results)
twas<-read.csv(paste0('/data/coxvgi/zhoud2/projects/gtex/asso/LDLq_',model,'_',tissue,'.csv'),header = T,stringsAsFactors = F)
twas<-twas[twas$pvalue<0.05,] #only run MR for TWAS-significant genes
gene_list<-intersect(gene_list,twas[,1])

#load ldsc
ldsc<-read.table('/data/coxvgi/zhoud2/data/gtex/geno/v8/ldsc/gtex_pruned.score.ld',header = T,stringsAsFactors = F)
ldsc<-ldsc[,c(1,8)]

#gene name annotation
gene_name_anno<-read.table('/data/coxvgi/zhoud2/anno/gencode/37/gencode.v32.GRCh37.txt',stringsAsFactors = F,header = T)
gene_name_anno<-gene_name_anno[,c('geneid','genename')]

#built output df
#output<-data.frame('geneid'=NA,'p_real'=NA,'beta_real'=NA,'lower_real'=NA,'upper_real'=NA,'p_shuffled'=NA,'beta_shuffled'=NA,'lower_shuffled'=NA,'upper_shuffled'=NA)
output<-data.frame('geneid'=NA)

##for test
#i=which(gene_list=='ENSG00000134243')
#i=which(gene_list=='ENSG00000198670')
#i=which(gene_list=='ENSG00000132170')
#i=which(gene_list=='ENSG00000107798')


#subjob start and end id
i_start=(args-1)*10+1
i_end=args*10
if (i_end>length(gene_list)){
  i_end=length(gene_list)
}


if (i_start<length(gene_list)){
  for (i in i_start:i_end){
    print(i)
    output[i,1]<-gene<-gene_list[i]
    gene_name<-gene_name_anno[which(gene_name_anno[,1]==gene_list[i]),2]
    
    #load genotype in dosage and allele info
    gtex_geno_path<-'/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/dosage_1m/'
    gtex_geno_info_path<-'/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/info_1m/'
    gtex_geno<-readRDS(paste0(gtex_geno_path,gene,'.rds'))
    gtex_info<-readRDS(paste0(gtex_geno_info_path,gene,'.rds'))
    
    #overlap snps with gwas results
    if(length(which(colnames(gtex_geno) %in% gwas$rsid))==0){next}
    gtex_geno<-gtex_geno[,c(1,which(colnames(gtex_geno) %in% gwas$rsid))]
    gtex_info<-gtex_info[which(gtex_info$rsid %in% gwas$rsid),c(1,4)]
    
    #load expression
    gtex_exp_path<-'/data/coxvgi/zhoud2/data/gtex/exp/v8/weighted/'
    gtex_exp<-readRDS(paste0(gtex_exp_path,tissue,'/',gene,'.rds'))
    
    #run eqtl, get beta_eqtl
    gtex_geno[,1]<-paste0('GTEX.',gtex_geno[,1])
    gtex_exp<-gtex_exp[gtex_exp$tissue==tissue,c(1,2)]
    d<-merge(gtex_exp,gtex_geno,by=1)
    eqtl<-data.frame('rsid'=gtex_info$rsid,'eqtl_beta'=NA,'eqtl_se'=NA,'eqtl_p'=NA,stringsAsFactors = F)
    for (k in 1:nrow(eqtl)){
      eqtl_fit<-summary(lm(d$exp~d[,k+2]))
      eqtl[k,2]<-eqtl_fit$coefficients[2,1]
      eqtl[k,3]<-eqtl_fit$coefficients[2,2]
      eqtl[k,4]<-eqtl_fit$coefficients[2,4]
    }
    eqtl<-merge(eqtl,gtex_info,by='rsid')
    
    #merge eqtl results with gwas results
    df<-merge(eqtl,gwas,by='rsid')
    
    #match eqtl and gwas effect alleles
    df$gwas_beta<-ifelse(df$counted_allele==df$effect_allele_gwas,df$beta_gwas,(df$beta_gwas*-1))
    
    #merge with ldsc
    df<-merge(df,ldsc,by=1)
    
    #skip genes with less than 20 obs (no way for cv)
    if (nrow(df)<20){next}
    
    #no of sig gwas loci (only consider these for pleiotropy control)
    n_gwas<-length(which(df$p_gwas<0.05))
    
    #identity matrix
    df<-df[order(df$p_gwas),]
    df<-df[,c('rsid','gwas_beta','eqtl_beta','ldscore')] 
    if(n_gwas>0){
      df[,(ncol(df)+1):(ncol(df)+n_gwas)]<-diag(nrow(df))[,1:n_gwas]
    }
    
    #---assign x and y---
    y=df[,'gwas_beta'];x=as.matrix(df[,3:ncol(df)])
    #scale
    y=scale(y)
    x=apply(x, 2, function(x) scale(x))
    
    #---real data bootstrap lasso---
    set.seed(i)
    fit<-try(RB_LASSO(x=x,y=y,intercept = T,standardize = T,nfolds = nf,alpha=0.05))
    if(!('try-error' %in% class(fit))){
      
      output[i,'p_real']<-fit$pvalue
      output[i,'beta_real']<-fit$beta[1]
      output[i,'lower_real']<-fit$interval[1,1]
      output[i,'upper_real']<-fit$interval[2,1]
    }
    
    
    #---data shuffling----
    set.seed(i+100)
    y<-sample(y,length(y),replace = F)
    
    #---shuffled data bootlasso find CI ----
    set.seed(i)
    fit<-try(RB_LASSO(x=x,y=y,intercept = T,standardize = T,nfolds = nf,alpha=0.05))
    if(!('try-error' %in% class(fit))){
      output[i,'p_shuffled']<-fit$pvalue
      output[i,'beta_shuffled']<-fit$beta[1]
      output[i,'lower_shuffled']<-fit$interval[1,1]
      output[i,'upper_shuffled']<-fit$interval[2,1]
    }
    
  }
  
  #sig or not according to 95%CI
  output$sig_real_ci<-ifelse(output$lower_real*output$upper_real>0,'sig','nonsig')
  output$sig_shuffled_ci<-ifelse(output$lower_shuffled*output$upper_shuffled>0,'sig','nonsig')
  
  #rm NA lines
  output<-output[!is.na(output[,2]),]
  
  #output
  write.csv(output,paste0('/data/coxvgi/zhoud2/projects/gtex/mr/rbl_noresample/p',args,'.csv'),quote = F,row.names = F)
  
}





















