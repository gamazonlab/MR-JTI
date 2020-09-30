args=as.numeric(commandArgs(TRUE))

#simulation times per gene
n_simulations=100
#n of causal genes
n_causal_genes=100
#sample size
n_samples=c(1000,5000,10000,25000,50000,100000,250000,500000)[args]
#heritability for all the causal genes
h2_all_gene=0.5

#generate true expression
set.seed(2020)
true_exp_matrix=matrix(rnorm(n_samples*n_causal_genes,0,1),ncol = 100)

#generate effect size of expression on trait
set.seed(2020)
beta<-rnorm(n_causal_genes,0,(h2_all_gene/n_causal_genes)^0.5)
#simulate phenotypes
set.seed(2021)
pheno<-true_exp_matrix %*% beta + rnorm(n_samples,mean = 0,sd=(1-h2_all_gene/n_causal_genes*n_causal_genes)^0.5)

#check
#ans<-summary(lm(pheno~true_exp_matrix %*% beta))


#load the real performance (r2 in exteral data) for each model 
perf<-read.table('/data/c***/z***/projects/gtex/simulation/power/summary.txt',header = T,stringsAsFactors = F)

#consider r2>0.04 genes as candidate causal genes.
perf$r_max<-apply(perf[,c('r_st','r_xt','r_ut')],1,function(x) max(x))
perf<-perf[perf$r_max>0.2,]

#sampling causal genes for simulation
set.seed(1)
perf<-perf[sample(seq(1,nrow(perf)),n_causal_genes,replace = F),]


#output df
output<-data.frame(sample_size=NA,gene_i=NA,simu_i=NA,p_st=NA,p_ut=NA,p_xt=NA)

i=1
for (gene_i in 1:n_causal_genes){
  print(gene_i)
  for (simu_i in 1:n_simulations){
    #predixcan
    set.seed(gene_i+simu_i+1)
    st_exp=rnorm(n_samples,0,1)
    #utmost
    set.seed(gene_i+simu_i+2)
    ut_exp=rnorm(n_samples,0,1)
    #jti
    set.seed(gene_i+simu_i+3)
    xt_exp=rnorm(n_samples,0,1)
    
    exp<-as.matrix(data.frame(true=true_exp_matrix[,gene_i],st=st_exp,ut=ut_exp,xt=xt_exp))
    
    #decorrelate
    c1 <- var(exp)  # find the current correlation matrix
    chol1 <- solve(chol(c1)) # cholesky decomposition to get independence
    exp_decor<-exp %*% chol1 
    
    #check independence
    #zapsmall(cor(exp_decor))
    
    perf_matrix <- matrix( 
      c(1, perf$r_st[gene_i], perf$r_ut[gene_i], perf$r_xt[gene_i], 
        perf$r_st[gene_i], 1, 0, 0,
        perf$r_ut[gene_i], 0, 1, 0,
        perf$r_xt[gene_i], 0, 0, 1), ncol=4 )
    
    chol2 <- try(chol(perf_matrix))
    if('try-error' %in% class(chol2)){next}
    
    #generate simulated expression levels
    exp_simu <- exp_decor %*% chol2 * sd(true_exp_matrix[,gene_i]) + mean(true_exp_matrix[,gene_i])
    colnames(exp_simu)<-c('true','st','ut','xt')
    
    output[i,'sample_size']<-n_samples
    output[i,'gene_i']<-gene_i
    output[i,'simu_i']<-simu_i
    output[i,'p_st']<-cor.test(pheno,exp_simu[,'st'])$p.value
    output[i,'p_ut']<-cor.test(pheno,exp_simu[,'ut'])$p.value
    output[i,'p_xt']<-cor.test(pheno,exp_simu[,'xt'])$p.value
    
    i=i+1
  }
}

write.table(output,paste0('/data/c***/z***/projects/gtex/simulation/power/result/',n_samples,'.txt'),quote = F,row.names = F,sep = '\t')

