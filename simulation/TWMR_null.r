#subjobs
args<-as.numeric(commandArgs(TRUE)) 

#load real GTEx eqtl results
eqtl<-read.table('/data/c***/z***/projects/gtex/mr/twmr/liver/eqtl.txt',header = T,stringsAsFactors = F)
eqtl<-eqtl[,c('gene_id','rsid','slope','slope_se','pval_nominal','ref_allele','eff_allele')]
colnames(eqtl)<-c('gene_id','rsid','beta_eqtl','se_eqtl','p_eqtl','ref_allele_eqtl','eff_allele_eqtl')

#load real TWAS results (only run MR for the list of genes with significant TWAS results)
twas<-read.csv(paste0('/data/c***/z***/projects/gtex/asso/LDLq_xt_Liver.csv'),header = T,stringsAsFactors = F)
#twas<-twas[twas$pvalue<0.05,]
#gene_list<-twas[,1] #1559 genes

#randomly selet 1000 genes
set.seed(2020)
gene_list<-sample(twas$gene,1000,replace = F)

#load snp position
bim<-read.table('/data/c***/z***/projects/gtex/mr/twmr/liver/unique_eqtl.bim',stringsAsFactors = F)
colnames(bim)<-c('chr','SNP','rm','pos','a1','a2')

#for test
#i=which(gene_list=='ENSG00000134243')
#i=which(gene_list=='ENSG00000163964')

#output df
output<-data.frame(gene_id=NA,alpha=NA,SE=NA,P=NA,Nsnps=NA,Ngene=NA)

#subjob start and end for each subjob
i_start=(args-1)*10+1
i_end=min((args)*10,length(gene_list))

for (i in i_start:i_end){
  print(i)
  gene_id<-gene_list[i]
  print(gene_id)
  
  #step 1. find independent eqtls to the focal gene 
  eqtl_unclumped<-eqtl[eqtl$gene_id==gene_id,]
  eqtl_unclumped<-eqtl_unclumped[,c('rsid','eff_allele_eqtl','ref_allele_eqtl','beta_eqtl','se_eqtl','p_eqtl')]
  colnames(eqtl_unclumped)<-c('SNP','eff_allele','ref_allele','BETA','se','P')
  eqtl_unclumped<-merge(eqtl_unclumped,bim,by='SNP')
  eqtl_unclumped<-eqtl_unclumped[,c('SNP','chr','pos','eff_allele','ref_allele','BETA','se','P')]
  #write the list of unclumped snps
  write.table(eqtl_unclumped,paste0('/data/c***/z***/projects/gtex/mr/twmr/liver/tmp/clumping/',gene_id),quote = F,row.names = F,sep='\t') 
  #LD clumping
  cmd<-paste0('plink --bfile /data/c***/z***/projects/gtex/mr/twmr/liver/unique_eqtl  --clump /data/c***/z***/projects/gtex/mr/twmr/liver/tmp/clumping/',gene_id,'  --clump-r2 0.1 --clump-field P --clump-p1 0.001 --clump-p2 0.001  --clump-kb 500  --out /data/c***/z***/projects/gtex/mr/twmr/liver/tmp/clumping/',gene_id) 
  system(cmd,wait = T)
  #load clumped SNPs
  if(!(file.exists(paste0('/data/c***/z***/projects/gtex/mr/twmr/liver/tmp/clumping/',gene_id,'.clumped')))){next}
  clumped_rsid<-read.table(paste0('/data/c***/z***/projects/gtex/mr/twmr/liver/tmp/clumping/',gene_id,'.clumped'),header = T,stringsAsFactors = F)  
  eqtl_1<-eqtl[eqtl$gene_id==gene_id,]
  eqtl_1<-eqtl_1[which(eqtl_1$rsid %in% clumped_rsid$SNP),]
  
  #step 2. find all other genes associated with the independent eQTLs to the focal gene
  genes_shared_eqtls<-eqtl[(eqtl$rsid %in% eqtl_1$rsid),'gene_id']
  
  #step 3. find all eQTLs of the gene list in step 2
  eqtl_2<-eqtl[which(eqtl$gene_id %in% genes_shared_eqtls),]
  #write tmp files for the included eqtls
  write.table(eqtl_2[,c(2,3)],paste0('/data/c***/z***/projects/gtex/mr/twmr/liver/tmp/info/',gene_id),quote = F,sep='\t',row.names = F)
  #extract these eqtls
  cmd=paste0('plink --bfile /data/c***/z***/projects/gtex/mr/twmr/liver/unique_eqtl --extract /data/c***/z***/projects/gtex/mr/twmr/liver/tmp/info/',gene_id,' --make-bed --out /data/c***/z***/projects/gtex/mr/twmr/liver/tmp/geno/',gene_id)
  system(cmd,wait = T,ignore.stdout=T,ignore.stderr=T)
  if(!(file.exists(paste0('/data/c***/z***/projects/gtex/mr/twmr/liver/tmp/geno/',gene_id,'.bed')))){next}
  #LD pruning
  cmd=paste0('plink --bfile /data/c***/z***/projects/gtex/mr/twmr/liver/tmp/geno/',gene_id,' --indep-pairwise 50 5 0.1 --out /data/c***/z***/projects/gtex/mr/twmr/liver/tmp/info/',gene_id)
  system(cmd,wait = T,ignore.stdout=T,ignore.stderr=T)
  #extract LD-pruned snps and transfer them to dosage data
  cmd=paste0('plink --bfile /data/c***/z***/projects/gtex/mr/twmr/liver/tmp/geno/',gene_id,' --extract /data/c***/z***/projects/gtex/mr/twmr/liver/tmp/info/',gene_id,'.prune.in --recode A --out /data/c***/z***/projects/gtex/mr/twmr/liver/tmp/pruned/',gene_id)
  system(cmd,wait = T,ignore.stdout=T,ignore.stderr=T)
  if(!(file.exists(paste0('/data/c***/z***/projects/gtex/mr/twmr/liver/tmp/pruned/',gene_id,'.raw')))){next}
  #load dosage data
  dosage<-read.table(paste0('/data/c***/z***/projects/gtex/mr/twmr/liver/tmp/pruned/',gene_id,'.raw'),header = T,stringsAsFactors = F)
  dosage_sample=dosage$FID
  dosage<-dosage[,-c(1,2,3,4,5,6)]
  
  #all the eqtls
  snp_list<-as.character(sapply(colnames(dosage),function(x) strsplit(x,"[_]")[[1]][1]))
  eqtl_2<-eqtl_2[which(eqtl_2$rsid %in% snp_list),]
  #all the eGenes
  genes_shared_eqtls<-unique(eqtl_2$gene_id)
  if(!(gene_id %in% genes_shared_eqtls)){next}
  
  #generate the matrix (input of twmr)
  mat<-as.data.frame(matrix(data=NA,ncol = length(genes_shared_eqtls)+1, nrow = length(snp_list)))
  colnames(mat)<-c('GENES',genes_shared_eqtls)
  mat[,'GENES']<-snp_list
  
  #simulate trait (samples used for eqtl analysis are filtered out, so eqtl and gwas are independent)
  set.seed(i)
  trait<-rnorm((nrow(dosage)-208),mean = 0,sd=1)
  
  #estimate beta for each snp-gene pairs and each snp-trait pairs
  for(j in 1:length(genes_shared_eqtls)){
    #load expression
    if(!(file.exists(paste0('/data/c***/z***/data/gtex/exp/v8/weighted/Liver/',genes_shared_eqtls[j],'.rds')))){next}
    exp<-readRDS(paste0('/data/c***/z***/data/gtex/exp/v8/weighted/Liver/',genes_shared_eqtls[j],'.rds'))
    exp<-exp[exp$tissue=='Liver',]
    exp$sampleid<-paste0('GTEX-',sapply(exp$sampleid,function(x) strsplit(x,"[.]")[[1]][2]))
    exp<-exp[,c('sampleid','exp')]
    #merge with dosage data
    dosage_tmp<-dosage
    dosage_tmp$sampleid=dosage_sample
    df<-merge(exp,dosage_tmp,by='sampleid')
    colnames(df)<-c('sampleid','exp',sapply(colnames(df)[-c(1,2)], function(x) strsplit(x,"[_]")[[1]][1]))
    
    for(k in 1:length(snp_list)){
      #estimate beta eqtl
      ans<-summary(lm(df$exp~df[,snp_list[k]]))
      if(ans$coefficients['df[, snp_list[k]]','Pr(>|t|)']<0.05){
        mat[k,genes_shared_eqtls[j]]<-ans$coefficients['df[, snp_list[k]]','Estimate']
      }
      #estimate beta gwas
      colnames(dosage_tmp)<-sapply(colnames(dosage_tmp), function(x) strsplit(x,"[_]")[[1]][1])
      ans<-summary(lm(trait~dosage_tmp[-which(dosage_tmp$sampleid %in% exp$sampleid),snp_list[k]]))
      mat[k,'BETA_GWAS']<-ans$coefficients[2,'Estimate']
    }
  }
  mat[is.na(mat)]<-0
  
  if(length(genes_shared_eqtls)>1){
    mat_target<-mat[,c('GENES',gene_list[i])]
    mat<-mat[,-c(1,which(colnames(mat)==gene_list[i]))]
    mat<-cbind(mat_target,mat)
  }
  
  
  #write matrix file
  write.table(mat,paste0('/data/c***/z***/projects/gtex/mr/twmr/liver/tmp/input/',gene_id,'.matrix'),quote = F,row.names = F,sep='\t')
  
  #generate ld cor matrix
  ld<-cor(dosage)
  #write ld cor matrix
  write.table(ld,paste0('/data/c***/z***/projects/gtex/mr/twmr/liver/tmp/input/',gene_id,'.ld'),quote = F,row.names = F,col.names = F,sep='\t')
  
  
  
  #----run TWMR-----
  
  main_dir='/data/c***/z***/projects/gtex/mr/twmr/liver/tmp/input'
  
  Ngwas<-nrow(dosage)-208  #208 is the sample size of eQTL of Liver, nrow(dosage) is the total number of samples for GTEx. Ngwas is the sample size of the simulated GWAS. So the GWAS samples and the eQTL samples are independent.
  N_eQTLs<-length(gene_list)
  out<-c("gene_id","alpha","SE","P","Nsnps","Ngene")
  
  file<-paste0(main_dir,'/',gene_id,".matrix")
  filecluster<-read.table(file,header=T)
  beta<-as.matrix(filecluster[,2:(length(filecluster[1,])-1)])
  
  x<-colSums(abs(beta))
  remove<-which(x==0)
  
  if(length(remove)>0) {
    beta<-beta[,-remove]}
  beta<-as.matrix(beta)
  
  gamma<-as.matrix(filecluster[,length(filecluster[1,])])
  
  LDmatrix<-paste0(main_dir,'/',gene_id,".ld")
  C<-read.table(LDmatrix,header=F)
  
  C<-as.matrix(C[,1:length(C[,1])])
  
  S<-try(t(beta)%*%solve(C)%*%beta)
  if('try-error' %in% class(S)){next}
  
  H<-try((1-1/sqrt(Ngwas))*S+(1/sqrt(Ngwas))*diag(length(S[,1])))  # we tested Ngwas = 3781 and also tested other GWAS sample sizes with consistent results
  if('try-error' %in% class(H)){next}
  
  alpha<-solve(H)%*%(t(beta)%*%solve(C)%*%gamma)
  
  alpha<-as.vector(alpha)
  
  C_inv <- solve(C)
  GCG_inv <- t(beta) %*% solve(C) %*% beta
  GCG_inv<-(1-1/sqrt(Ngwas))*GCG_inv+(1/sqrt(Ngwas))*diag(length(GCG_inv[,1])) # we tested Ngwas = 3781 and also tested other GWAS sample sizes with consistent results
  GCG_inv<-solve(GCG_inv)
  
  df_dg <- GCG_inv %*% t(beta) %*% C_inv
  df_dG <- (GCG_inv %x% (t(gamma) %*% C_inv %*% ((beta %*% GCG_inv %*% t(beta)) %*% C_inv + diag(nrow(beta))))) + ((-t(gamma) %*% C_inv %*% beta %*% GCG_inv) %x% (GCG_inv %*% t(beta) %*% C_inv))
  J <- cbind(df_dG, df_dg)
  
  SEs<-c(rep(1/sqrt(length(gene_list)),length(beta[1,])*length(beta[,1])),rep(1/sqrt(Ngwas),length(gamma[,1])))
  R<-diag(length(beta[1,])+1)
  Sigma <- (SEs %*% t(SEs)) * (C %x% R)   
  V <- J %*% Sigma %*% t(J)
  se<- sqrt(V[1,1])
  
  N=length(beta[,1])
  Ngene=length(beta[1,])
  Z<-alpha[1]/se
  pval<-2*pnorm(abs(Z),lower.tail=FALSE)
  line<-c(colnames(filecluster)[2],alpha[1],se,pval,N,Ngene)
  out<-rbind(out,line)
  
  output[i,1:6]<-line
  
  
}

output<-output[!is.na(output[,1]),]

write.table(output,file=paste0('/data/c***/z***/projects/gtex/mr/twmr/liver/t1e/0.05_',args,'.txt'),quote=F,sep='\t',row.names=F)













