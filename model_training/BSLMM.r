args<-as.numeric(commandArgs(TRUE))

#set up sub jobs
run_id<-1
run_list<-list()
for (i in c(13,22)){  #tissues 13:bfc and 22:lcl
  for (j in 1:1400){  #subjobs per tissue
    run_list[[run_id]]<-c(i,j)
    run_id=run_id+1
  }
}
run_i<-run_list[[args]]

#tissue list
tissue_list<-dir('/data/c***/z***/data/gtex/exp/v8/weighted/')
tissue<-tissue_list[run_i[1]]

#make dirs
dir.create(paste0('/data/c***/z***/projects/gtex/weights/bslmm/',tissue))
dir.create(paste0('/data/c***/z***/projects/gtex/bslmm/fusion_output/',tissue))
dir.create(paste0('/data/c***/z***/projects/gtex/bslmm/tmp/args/',tissue,'_',run_i[2]))

#gene list
gene_list<-dir(paste0('/data/c***/z***/data/gtex/exp/v8/weighted/',tissue))
genotype_list<-dir('/data/c***/z***/data/gtex/geno/v8/gene/dosage_1m/')
gene_list<-intersect(gene_list,genotype_list)

#annotation for genes
anno<-read.table('/data/c***/z***/anno/gencode/37/gencode.v32.GRCh37.txt',header=T,stringsAsFactors = F)
anno<-anno[(anno$genetype %in% c('protein_coding','lncRNA','miRNA')),]
anno<-anno[(anno$chr %in% paste0('chr',seq(1,22))),]

#subjob start and end id
i_start=(run_i[2]-1)*15+1
i_end=run_i[2]*15
if (i_end>length(gene_list)){
  i_end=length(gene_list)
}

if (i_start<length(gene_list)){
  for (i in i_start:i_end){
    #print(i)
    
    geneid=sub('....$','',gene_list[i])
    
    #geneid='ENSG00000134243'
    #geneid='ENSG00000164308'
    
    #check exist
    if(file.exists(paste0('/data/c***/z***/projects/gtex/weights/bslmm/',tissue,'/',geneid,'.rds'))){next}
    
    #find chr pos
    chr=sub('^...','',anno[which(anno$geneid==geneid),'chr'])
    left=max(anno[which(anno$geneid==geneid),'left']-500000,1)
    right=anno[which(anno$geneid==geneid),'right']+500000
    
    #extract genotype
    cmd<-paste0('plink --bfile /data/c***/z***/data/gtex/geno/v8/liftover/pruning/tissue_specific/',tissue,' --chr ',chr,' --from-bp ',left,' --to-bp ',right,' --make-bed --out /data/c***/z***/projects/gtex/bslmm/tmp/geno/',tissue,'_',geneid)
    system(cmd,wait = T,ignore.stdout=T,ignore.stderr=T)
    
    #load expression data
    exp<-readRDS(paste0('/data/c***/z***/data/gtex/exp/v8/weighted/',tissue,'/',geneid,'.rds'))
    exp<-exp[which(exp$tissue==tissue),]
    exp$sampleid<-sub('GTEX.','GTEX-',exp$sampleid)
    
    #load fam file 
    fam<-read.table(paste0('/data/c***/z***/projects/gtex/bslmm/tmp/geno/',tissue,'_',geneid,'.fam'),stringsAsFactors = F)
    #generate the trait (expression, here) file by modfify the fam file
    for (j in 1:nrow(fam)){
      if(fam[j,'V1'] %in% exp$sampleid){
        fam[j,'V6']<-exp[which(exp$sampleid==fam[j,'V1']),'exp']
      }
    }
    
    #write fam file as the trait file (will be used for bslmm)
    write.table(fam,paste0('/data/c***/z***/projects/gtex/bslmm/tmp/geno/',tissue,'_',geneid,'.fam'),quote = F,row.names = F,col.names = F,sep = '\t')
    
    #set dir
    setwd(paste0('/data/c***/z***/projects/gtex/bslmm/tmp/args/',tissue,'_',run_i[2]))
    
    #creat a link (gemma write output here)
    if(!file.exists(paste0('/data/c***/z***/projects/gtex/bslmm/tmp/args/',tissue,'_',run_i[2],'/output'))){
      cmd=paste0(paste0('ln -s ./ output'))
      system(cmd,wait = T)
    }
    
    #run 5x BSLMM
    cmd=paste0('Rscript /data/c***/z***/tools/fusion/fusion_twas-master/FUSION.compute_weights.R  --bfile /data/c***/z***/projects/gtex/bslmm/tmp/geno/',tissue,'_',geneid,'  --tmp ',geneid,'  --out /data/c***/z***/projects/gtex/bslmm/fusion_output/',tissue,'/',geneid,'  --models top1,bslmm,enet --PATH_gemma /data/c***/z***/tools/bslmm/gemma-0.98.1 --PATH_gcta /home/zhoud2/tools/gcta/gcta64 --hsq_p 0.1')
    system(cmd,wait = T)
    
    #clean plink files
    cmd=paste0('rm /data/c***/z***/projects/gtex/bslmm/tmp/geno/',tissue,'_',geneid,'*')
    system(cmd,wait = T)
    
    #load results
    if(!file.exists(paste0('/data/c***/z***/projects/gtex/bslmm/fusion_output/',tissue,'/',geneid,'.wgt.RDat'))){next}
    load(paste0('/data/c***/z***/projects/gtex/bslmm/fusion_output/',tissue,'/',geneid,'.wgt.RDat'))
    snps$weight=wgt.matrix[,'bslmm']
    snps$r2=r2=cv.performance['rsq','bslmm']
    snps$p=p=cv.performance['pval','bslmm']
    snps$gene=geneid
    colnames(snps)[1:6]<-c('chr','rsid','drop','pos','counted_allele','ref_allele')
    snps$chr_pos<-paste0(snps$chr,'_',snps$pos)
    snps<-snps[,c('gene','rsid','chr_pos','ref_allele','counted_allele','weight','r2','p')]
    snps<-snps[which(snps$weight!=0),]
    
    #write results
    if(nrow(snps)>0 & r2>0.01 & p<0.05){
      saveRDS(snps,paste0('/data/c***/z***/projects/gtex/weights/bslmm/',tissue,'/',geneid,'.rds'))
    }
    
    
  }
}









