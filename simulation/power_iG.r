
args<-commandArgs(TRUE) #[1] tissue id [2] model

info<-read.table('/data/coxvgi/zhoud2/projects/cross_tissue/info/gtex_info.txt',header = T,stringsAsFactors = F)
tissue_list<-info[,1]
tissue=tissue_list[as.numeric(args[1])]
model=args[2]  #xt and st

#genelist
gene_list_st<-dir(paste0('/data/coxvgi/zhoud2/data/biovu/pred_exp_23k/st/',tissue))
gene_list_xt<-dir(paste0('/data/coxvgi/zhoud2/data/biovu/pred_exp_23k/xt/',tissue))
gene_list_all<-intersect(gene_list_xt,gene_list_st)

#n_iG
n_iG<-length(dir(paste0('/data/coxvgi/zhoud2/projects/gtex/weights/',model,'/',tissue)))

#filter samples without ldl
ldl<-read.table('/data/coxvgi/zhoud2/data/biovu/pheno/lipid/multixcan_pheno_ldl.txt',header = T,stringsAsFactors = F)
ldl<-ldl[!is.na(ldl[,2]),]
sd_ldl<-sd(ldl[,2])

#100 random genes
set.seed(args[1])
gene_list<-sample(gene_list_all,100,replace = F)

out<-as.data.frame(matrix(data=NA,ncol=4,nrow=0))
colnames(out)<-c('gene','X0.10','X0.30','X0.50')

for (i in 1:length(gene_list)){
  print(i)
  gene=gene_list[i]
  d_weight<-readRDS(paste0('/data/coxvgi/zhoud2/projects/gtex/weights/',model,'/',tissue,'/',gene))
  
  d_exp<-readRDS(paste0('/data/coxvgi/zhoud2/data/biovu/pred_exp_23k/',model,'/',tissue,'/',gene))
  X_i=1
  for (X in c(0.10,0.30,0.50)){
    #X=0.05
    err_sd<-1-X^2*var(d_exp[,2])/sd_ldl^2
    
    p<-c()
    for (k in 1:500){
      #print(k)
      d_exp$y<-X*d_exp$pred+rnorm(nrow(d_exp),mean=0,sd=err_sd)
      fit<-summary(lm(d_exp$y~d_exp$pred))
      p[k]<-fit$coefficients[2,4]
    }
    out[i,1]<-gene
    out[i,X_i+1]<-length(which(p<0.05/n_iG))/500
    X_i=X_i+1
  }
}

write.table(out,paste0('/data/coxvgi/zhoud2/projects/gtex/simulation/power_iG/',args[2],'/',tissue,'.txt'),quote = F,row.names = F,sep='\t')

