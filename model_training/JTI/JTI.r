
#-----Joint-Tissue-gene-expression-Imputation-(JTI)-----
cat('\n---Joint-Tissue-gene-expression-Imputation-(JTI)---\n')

library('glmnet')
library("optparse")

option_list = list(
  make_option("--tissue", action="store", default=NA, type='character',
              help="tissue name [required]"),
  make_option("--geneid", action="store", default=NA, type='character',
              help="ENSG geneid [required]"),
  make_option("--plink_path", action="store", default="plink", type='character',
              help="Path to plink [%default]"),
  make_option("--genotype_path", action="store", default=NA, type='character',
              help="genotype file in plink bfile format '.bed, .bim, .fam' [required]"),
  make_option("--expression_path", action="store", default=NA, type='character',
              help="path to expression file [required]"),
  make_option("--tmp_folder", action="store", default=NA, type='character',
              help="tmp folder for raw dosage data [required]"),
  make_option("--gencode_path", action="store", default=NA, type='character',
              help="path to gencode annotation file [required]"),
  make_option("--out_path", action="store", default=NA, type='character',
              help="folder for output [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

#input
tissue=opt$tissue
geneid=opt$geneid
plink_path=opt$plink_path
genotype_path=opt$genotype_path
tmp_folder=opt$tmp_folder
expression_path=opt$expression_path
gencode_path=opt$gencode_path
out_path=opt$out_path

cat(' INFO loading gene position annotation ...\n')
#load gene annotation file
gencode<-read.table(gencode_path,header = T,stringsAsFactors = F)

cat(' INFO mkdir tmp folders ...\n')
#mkdir tmp folder. will be cleaned
options(warn=-1)
dir.create(tmp_folder)
dir.create(paste0(tmp_folder,'/',tissue,'_',geneid))
tmp_folder=paste0(tmp_folder,'/',tissue,'_',geneid)
options(warn=0)

#get chr pos for the gene
chr=as.numeric(sub('^...','',gencode[which(gencode$geneid==geneid),'chr']))
pos_from=max(1,gencode[which(gencode$geneid==geneid),'left']-1000000)
pos_to=gencode[which(gencode$geneid==geneid),'left']+1000000

cat(' INFO generating dosage genotype data ...\n')
#extract genotypes from plink file to dosage file
cmd=paste0(plink_path,' --bfile ',genotype_path,' --chr ',chr,' --from-bp ',pos_from,' --to-bp ',pos_to,' --recode A --out ',tmp_folder,'/',geneid)
system(cmd,ignore.stdout=T,ignore.stderr=T)
cmd=paste0(plink_path,' --bfile ',genotype_path,' --chr ',chr,' --from-bp ',pos_from,' --to-bp ',pos_to,' --make-bed --out ',tmp_folder,'/',geneid)
system(cmd,ignore.stdout=T,ignore.stderr=T)

#load dosage file
dosage_raw<-try(read.table(paste0(tmp_folder,'/',geneid,'.raw'),header = T,stringsAsFactors = F))
if('try-error' %in% class(dosage_raw)){stop('no SNP available for this gene')}
dosage<-dosage_raw[,-c(1,3:6)]  #rm some cols
colnames(dosage)<-c('sampleid',sapply(colnames(dosage)[-1],function(x) strsplit(x,"[_]")[[1]][1])) #rm the counted allele from rsnum
dosage[,-1]<-round(apply(dosage[,-1], 2, function(x) ifelse(is.na(x),mean(x,na.rm=T),x)),3) #post imputation imputation. NA replaced by mean

#load SNP info (alleles)
snp_info_raw<-read.table(paste0(tmp_folder,'/',geneid,'.bim'),stringsAsFactors = F)
snp_info_raw$counted_allele<-snp_info_raw$V5
snp_info_raw$ref_allele<-snp_info_raw$V6
snp_info_raw$chr_bp<-paste0(snp_info_raw$V1,'_',snp_info_raw$V4)
colnames(snp_info_raw)[2]<-'rsid'
snp_info<-snp_info_raw[,c('rsid','chr_bp','ref_allele','counted_allele')]

#load expression
cat(' INFO loading expression data ...\n')
exp<-read.table(expression_path,header = T,stringsAsFactors = F)
    

#model training
cat(' INFO initializing model training ...\n')
#fit single tissue model to get lambda range 
exp_st<-exp[exp$tissue==tissue,] #load single tissue data
d_st<-merge(exp_st,dosage,by='sampleid') #merge
set.seed(as.numeric(sub('^....','',geneid))) 
fit_st<-cv.glmnet(x=as.matrix(d_st[,6:ncol(d_st)]),y=as.matrix(d_st[,'exp']), nfolds = 5,keep = T,alpha=0.5) 
lambda_list<-fit_st$lambda[max(1,which(fit_st$lambda==fit_st$lambda.min)-10):min(length(fit_st$lambda),which(fit_st$lambda==fit_st$lambda.min)+10)] 

#generate dataframe for grid search and performance collection
r_matrix<-as.data.frame(matrix(data=NA,nrow=(4*4),ncol=6))
i_loop<-1
for (exp_power in c(1,4,16,64)){
  for (dhs_power in c(1,4,16,64)){
    r_matrix[i_loop,1:2]<-c(exp_power,dhs_power)
    i_loop=i_loop+1
  }
}
colnames(r_matrix)<-c('exp_power','dhs_power','lambda','window','r_test','p_test')
r_matrix$window='1m'

#sample id map for each fold (fold=5)
d_tmp<-merge(exp,dosage,by='sampleid')
sample_all<-unique(d_tmp$sampleid)
set.seed(as.numeric(sub('^....','',geneid))+1)
id_map<-data.frame(sampleid=sample_all,id_map=sample(rep(seq(1,5),ceiling(length(sample_all)/5)),length(sample_all)),stringsAsFactors = F)

#beta list of each hyper-parameter pairs
beta_list<-list()

cat(' INFO grid searching ...\n')
#grid search for hyper-parameter pairs
for (j in 1:nrow(r_matrix)){
  exp_power=r_matrix[j,1];dhs_power=r_matrix[j,2]
  exp$w<-(exp$exp_w)^exp_power*(exp$dhs_w)^dhs_power
  
  #merge
  d<-merge(exp,dosage,by='sampleid')
  
  #rm tissues with low similarity levels (w<0.1)
  d<-d[d[,'w']>=0.1,]; d<-d[!is.na(d[,1]),]
  
  #id map
  d<-merge(id_map,d,by='sampleid')
  
  #target tissue position (the performance will be only estimated in the target tissue)
  tt_pos<-which(d$tissue==tissue)
  
  #cross-tissue weighted elastic net
  set.seed(as.numeric(sub('^....','',geneid))+2)
  ans<-try(cv.glmnet(x=as.matrix(d[,8:ncol(d)]),y=as.matrix(d[,'exp']),weights=d[,'w'],foldid=d[,'id_map'],lambda=lambda_list,keep = T,pmax=200)) 
  if ('try-error' %in% class(ans)){
    r_matrix[j,'r_test']<-0
  }else{
    #correlation between pred and obs expression levels
    cor_ans<-cor.test(ans$fit.preval[tt_pos,which(ans$lambda==ans$lambda.min)],d[tt_pos,'exp'])
    
    r_matrix[j,'r_test']<-cor_ans$estimate
    r_matrix[j,'p_test']<-cor_ans$p.value
    r_matrix[j,'lambda']<-ans$lambda.min
    
    #collect the weights
    beta_list[[j]]<-as.numeric(ans$glmnet.fit$beta[,which(ans$lambda==ans$lambda.min)])
    
  }
}

#find the best hyperparameters with the largest r_test
best_row=which.max(r_matrix[,'r_test'])[1]
r_test=r_matrix[best_row,'r_test']
p_test=r_matrix[best_row,'p_test']

cat(' INFO the r = ',r_test,', p = ',p_test,' \n')

#---output---
if (p_test<0.05 & r_test>0.1){
  snp_info$gene<-geneid
  snp_info$r2<-r_test^2
  
  snp_info$p<-p_test
  snp_info$lambda<-r_matrix[best_row,'lambda']
  snp_info$weight<-beta_list[[best_row]]
  snp_info<-snp_info[,c(5,1:4,9,6:8)]
  weight_file<-snp_info[snp_info$weight!=0,]
  if (nrow(weight_file)>0){
    write.table(weight_file,paste0(out_path,'/',geneid,'_',tissue,'.txt'),quote = F,row.names = F,sep = '\t')
  }
}

#cleaning
cat(' INFO cleaning tmp folder \n')
cmd=paste0('rm -r ',tmp_folder)
system(cmd,wait = T)

cat(' INFO done \n')



