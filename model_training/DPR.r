#DPR environment 
#ml GCC/6.4.0-2.28 Python/3.6.3 OpenMPI/2.1.1 numpy/1.13.1-Python-3.6.3 pandas/0.18.1-Python-3.6.3 scipy/0.19.1-Python-3.6.3  tabix/0.2.6 R
#source /data/c***/z***/tools/TIGAR/tiger_env/bin/activate
#export PYTHONPATH=$(python -c 'import sys; print(sys.path[-1])'):${PYTHONPATH}

#---R script for DPR model training---
args<-as.numeric(commandArgs(TRUE))

#subjob id
subjob=args

#training tissue
tissue='Brain_Frontal_Cortex_BA9' 

#load annotation for gene position
anno<-read.table('/data/c***/z***/anno/gencode/37/gencode.v32.GRCh37.txt',header = T,stringsAsFactors = F)
anno<-anno[,c('chr','left','right','geneid','genename')]

#load expression (residual)
exp_all<-read.table(paste0('/data/c***/z***/data/gtex/exp/v8/exp_residual/',tissue,'.v8.exp.residual'),header = T,stringsAsFactors = F)
colnames(exp_all)<-sapply(colnames(exp_all), function(x) sub('GTEX.','GTEX-',x))
exp<-merge(anno,exp_all,by='geneid')
exp$chr=as.numeric(sub('chr','',exp$chr))
exp<-exp[order(exp$chr),]

#assign 50 genes for each subjob (the 50 genes have to be on the same chromosome)
sub_id=0
for(i in 1:22){
  chr_tmp<-exp[exp$chr==i,1:4]
  chr_times<-ceiling(nrow(chr_tmp)/50)
  chr_tmp$subjob=rep((sub_id+1):(sub_id+chr_times),50)[1:nrow(chr_tmp)]
  sub_id=sub_id+chr_times
  if(i==1){
    chr_df=chr_tmp
  }else{
    chr_df<-rbind(chr_df,chr_tmp)  
  }
}
#check DO NOT RUN.
#chr_df<-chr_df[order(chr_df$subjob),]

exp<-exp[which(chr_df$subjob==subjob),]
chr=exp[1,'chr']

#mkdir
dir.create(paste0('/data/g***/z***/gtex/model_training/DPR/raw/',tissue,'/',subjob))
setwd(paste0('/data/g***/z***/gtex/model_training/DPR/raw/',tissue,'/',subjob))
dir.create('./output')

#write expression file (input of DPR)
exp<-exp[,c(2,3,4,1,5,seq(6,ncol(exp)))]
colnames(exp)[1:5]<-c('CHROM','GeneStart','GeneEnd','TargetID','GeneName')
write.table(exp,paste0('./exp.txt'),quote = F,sep='\t',row.names = F)
#write sample file (input of DPR)
sampleid_df<-as.data.frame(colnames(exp)[6:ncol(exp)],stringsAsFactors = F)
write.table(sampleid_df,paste0('./sampleid.txt'),row.names = F,col.names = F,sep='\t',quote = F)

### Train DPR model
data_dir=paste0('/data/g***/z***/gtex/model_training/DPR/raw/',tissue,'/',subjob)

setwd('/data/c***/z***/tools/TIGAR/TIGAR-master/')

cmd=paste0('/data/c***/z***/tools/TIGAR/TIGAR-master/TIGAR_Model_Train.sh --model DPR --Gene_Exp ',data_dir,'/exp.txt --train_sampleID ',data_dir,'/sampleid.txt --genofile /data/c***/z***/data/gtex/geno/v8/liftover/pruning/vcf/gtex_rsid_prune0.9.vcf.gz --chr ',chr,' --genofile_type vcf --Format GT  --maf 0.05 --hwe 0.00001 --cvR2 1 --dpr 1 --ES fixed --out_dir ',data_dir,'/output')
system(cmd,wait = T,ignore.stdout=T,ignore.stderr=T)


