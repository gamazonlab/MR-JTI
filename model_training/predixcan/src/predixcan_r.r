
cat('-----PrediXcan-R----- \n-----2021-07-27---- \n')

library(optparse)
library(glmnet)
library(RSQLite)
library(data.table)  #new added
library(foreach)   #new added
library(doParallel)   #new added

option_list = list(
  #general
  make_option("--main_dir", action="store", default=NA, type='character'),
  make_option("--plink_file_name", action="store", default=NA, type='character'),
  make_option("--expression_file_name", action="store", default=NA, type='character'),
  make_option("--annotation_file_name", action="store", default=NA, type='character'),
  make_option("--parallel", action="store_true", default=FALSE,
              help="max(n-1,1) core(s) will be used for parallel model training / association test, where n is the number of available cores"),
  
  #model training
  make_option("--model_training", action="store_true", default=FALSE),
  make_option("--subjob_id", action="store", default=1, type='integer'),
  make_option("--n_folds", action="store", default=5, type='integer',
              help="Number of cross-validation folds"),
  make_option("--n_genes_for_each_subjob", action="store", default=200, type='integer',
              help="Number of cross-validation folds"),
  make_option("--cis_window_size", action="store", default=1e6, type='integer'),
  make_option("--include_all_gene_types", action="store_true", default=FALSE,
              help = "otherwise, only protein_coding lincRNA, and miRNA will be included"),
  make_option("--only_include_variant_with_a_rsid", action="store_true", default=FALSE,
              help = "only variants with a name start with 'rs' will be included"),
  
  #generate db and cov
  make_option("--generate_db_and_cov", action="store_true", default=FALSE),
  make_option("--output_file_name", action="store", default=NA, type='character'),
  make_option("--clean_tmp", action="store_true", default=FALSE),
  
  
  #association test
  make_option("--asso_test", action="store_true", default=FALSE),
  make_option("--db_path", action="store", default=NA, type='character',
              help="Path to the weight file (.db) [required]"),
  make_option("--cov_path", action="store", default=NA, type='character',
              help="Path to the covariance file [required]"),
  make_option("--gwas_path", action="store", default=NA, type='character',
              help="Path to the GWAS summary statistics [required]"),
  make_option("--gene_list_path_for_asso_test", action="store", default=NA, type='character',
              help="Path to a file to specify a list of geneids need to run."),
  make_option("--asso_out_path", action="store", default=NA, type='character',
              help="Path to result file [required]"),
  
  make_option("--gwas_variant_col", action="store", default=NA, type='character',
              help="column name of gwas variant id [required]"),
  make_option("--gwas_beta_col", action="store", default=NA, type='character',
              help="column name of gwas beta [required, if OR is not provided]"),
  make_option("--gwas_or_col", action="store", default=NA, type='character',
              help="column name of gwas OR [required, if beta is not provided]"),
  make_option("--gwas_se_col", action="store", default=NA, type='character',
              help="column name of gwas se(beta) [one the three (se, pvalue, and zscore) is required]"),
  make_option("--gwas_zscore_col", action="store", default=NA, type='character',
              help="column name of gwas zscore [one the three (se, pvalue, and zscore) is required]"),
  make_option("--gwas_p_col", action="store", default=NA, type='character',
              help="column name of gwas pval [one the three (se, pvalue, and zscore) is required]"),
  make_option("--gwas_eff_allele_col", action="store", default=NA, type='character',
              help="column name of gwas effect allele [required]"),
  make_option("--gwas_ref_allele_col", action="store", default=NA, type='character',
              help="column name of gwas reference allele [required]")
  

)


opt = parse_args(OptionParser(option_list=option_list))


if(opt$model_training+opt$generate_db_and_cov+opt$asso_test != 1){
  stop('Please choose one of the three options "--model_training", "--generate_db_and_cov", or "--asso_test"')
}

#functions

model_training<-function(main_dir,
                         plink_file_name,
                         expression_file_name,
                         annotation_file_name,
                         cis_window_size=1e6,
                         i,
                         n_genes_for_each_subjob){
  
  #main_dir<-'/data/coxvgi/zhoud2/tools/predixcan'
  #cis_window_size=1e6
  #plink_file_name<-'geuvadis421'
  #expression_file_name<-'geuvadis_residual.txt'
  #annotation_file_name<-'gencode.v32.GRCh37.txt'
  #i=1
  
  #mkdir
  if(!dir.exists(paste0(main_dir,'/gene/'))){
    dir.create(paste0(main_dir,'/gene/'))
  }
  if(!dir.exists(paste0(main_dir,'/weights/'))){
    dir.create(paste0(main_dir,'/weights/'))
  }
  if(!dir.exists(paste0(main_dir,'/cov/'))){
    dir.create(paste0(main_dir,'/cov/'))
  }
  
  #load annotation
  print('INFO loading annotation ...')
  anno<-read.table(paste0(annotation_file_name),header = T,stringsAsFactors = F) #keeps the colnames consistant with the example file
  
  if(!opt$include_all_gene_types){
    anno<-anno[(anno$genetype %in% c('protein_coding','lincRNA','miRNA')),]
  }
  anno<-anno[(anno$chr %in% paste0('chr',seq(1,22))),]
  
  #load expression data
  print('INFO loading expression data ...')
  exp_all<-read.table(paste0(expression_file_name),header = T,stringsAsFactors = F)  #keeps IID in genotype file match with IID in expression data
  #gene_list
  gene_list<-intersect(exp_all$geneid,anno$geneid)
  
  i_start<-(i-1)*n_genes_for_each_subjob+1  
  i_end<-min(i*n_genes_for_each_subjob,length(gene_list))
  

  if (i_start<=length(gene_list)){
    print(paste0('INFO total gene no ',length(gene_list),', for this subjob, from ',i_start,' to ',i_end,' ...'))
    
    
    if(opt$parallel){
      
      registerDoParallel(cores = max(detectCores()-1, 1))

      foreach(i = i_start:i_end) %dopar% {
        
        geneid=gene_list[i] #geneid='ENSG00000003402'
        print(paste0('INFO processing gene no.',i,' ',geneid))
        
        #get gene position  #use same build
        chr=sub('^...','',anno$chr[which(anno$geneid==geneid)])
        pos_l=max(anno$left[which(anno$geneid==geneid)]-cis_window_size,1)
        pos_r=anno$right[which(anno$geneid==geneid)]+cis_window_size
        
        #generate expression df
        exp<-data.frame('IID'=as.character(colnames(exp_all)[-1]),'exp'=as.numeric(exp_all[exp_pos<-which(exp_all[,1]==geneid),][-1]))
        
        each_gene(main_dir=main_dir,plink_file_name=plink_file_name,chr=chr,pos_l=pos_l,pos_r=pos_r,geneid=geneid,exp=exp)

      }

      stopImplicitCluster()
      
    }else{
      
      for (i in i_start:i_end){
        geneid=gene_list[i] #geneid='ENSG00000003402'
        print(paste0('INFO processing gene no.',i,' ',geneid))
        
        #get gene position  #use same build
        chr=sub('^...','',anno$chr[which(anno$geneid==geneid)])
        pos_l=max(anno$left[which(anno$geneid==geneid)]-cis_window_size,1)
        pos_r=anno$right[which(anno$geneid==geneid)]+cis_window_size
        
        #generate expression df
        exp<-data.frame('IID'=as.character(colnames(exp_all)[-1]),'exp'=as.numeric(exp_all[exp_pos<-which(exp_all[,1]==geneid),][-1]))
        
        each_gene(main_dir=main_dir,plink_file_name=plink_file_name,chr=chr,pos_l=pos_l,pos_r=pos_r,geneid=geneid,exp=exp)
        
      }
      
    }

  }else{
    print('start_i > N of genes')
  }
  
}

  
  
snpid_helper = function(str){
  elements = strsplit(str,"[_]")[[1]]
  snpid = list()
  snpid$id = paste(elements[-length(elements)],collapse = '_')
  snpid$allele = elements[length(elements)]
  return(snpid)
}
  


each_gene<-function(main_dir,plink_file_name,chr,pos_l,pos_r,geneid,exp){
  #generate genotype dosage file for each gene
  #extract from plink
  cmd<-paste0('plink --bfile ',plink_file_name,' --chr ',chr,' --from-bp ',pos_l,' --to-bp ',pos_r,' --recodeA --out ',main_dir,'/gene/',geneid)
  system(cmd,ignore.stdout=T,ignore.stderr=T,wait=T)
  
  #load dosage
  dosage_raw<-try(read.table(paste0(main_dir,'/gene/',geneid,'.raw'),header = T,stringsAsFactors = F))
  
  if(opt$only_include_variant_with_a_rsid){
    dosage_raw<-dosage_raw[,c(seq(1,6),which(substr(colnames(dosage_raw),1,2)=='rs'))]
  }

  if(!('try-error' %in% class(dosage_raw))){
    #rm raw file
    cmd<-paste0('rm ',main_dir,'/gene/',geneid,'*'); system(cmd,wait = T)
    #del useless cols
    dosage<-dosage_raw[,-c(1,3:6)] 
    
    #check if there is a '-' to '.' issue in colnames
    if(length(intersect(exp$IID, dosage$IID))==0){
      #if so, replace all the . with -
      exp$IID = gsub('\\.','-',exp$IID)
    }
    
    #merge genotype and expression
    df<-merge(exp,dosage,by='IID')
    
    #run elastic net
    y<-as.matrix(df[,2])
    x<-as.matrix(df[,c(3:ncol(df))])
    rm(df)
    
    #skip genes with only one SNP
    if(ncol(x)>=2){
      
      set.seed(as.numeric(sub('^....','',geneid)))
      fit<-try(cv.glmnet(x=x,y=y, nfolds = 5,keep = T,alpha=0.5,nlambda=50,pmax=200))
      if(!('try-error' %in% class(fit))){
        beta=as.numeric(fit$glmnet.fit$beta[,which(fit$lambda==fit$lambda.min)])
        
        performance<-cor.test(y,fit$fit.preval[,which(fit$lambda == fit$lambda.min)])
        r<-as.numeric(performance$estimate)
        p<-as.numeric(performance$p.value)
        
        #generate weight df
        weights_df<-data.frame('rsid'=colnames(dosage)[-1],'weights'=beta,stringsAsFactors = F,'r2'=r^2,'p'=p)
        weights_df$effect_allele<-sapply(weights_df$rsid, function(x) snpid_helper(x)[['allele']])
        weights_df$rsid<-sapply(weights_df$rsid, function(x) snpid_helper(x)[['id']])
        
        #rm weight=0 rows
        weights_df<-weights_df[weights_df$weights!=0,]
        
        if (nrow(weights_df)>0 & r>0.1 & p<0.05){
          
          #generate covariance matrix
          colnames(dosage)[-1]<-sapply(colnames(dosage)[-1], function(x) snpid_helper(x)[['id']])
          dosage<-dosage[,c(which(colnames(dosage) %in% weights_df$rsid))]
          cov_df<-as.data.frame(matrix(data=NA,ncol=4,nrow=1)) 
          colnames(cov_df)<-c('GENE','RSID1','RSID2','VALUE')
          
          #calculate covariance
          if(nrow(weights_df)>1){
            snp_list<-colnames(dosage)
            o_i=1
            for (k in 1:length(snp_list)){
              for (j in k:length(snp_list)){
                cov_df[o_i,2]<-snp_list[k]
                cov_df[o_i,3]<-snp_list[j]
                cov_df[o_i,4]<-round(cov(dosage[,which(colnames(dosage)==snp_list[k])],dosage[,which(colnames(dosage)==snp_list[j])]),3)
                o_i=o_i+1
              }
            }
          }else{
            cov_df[,'RSID1']=cov_df[,'RSID2']=weights_df[1,'rsid']
            cov_df[,'VALUE']=cov(dosage,dosage)
          }
          
          cov_df[,1]<-geneid
          
          #output cov
          write.table(cov_df,paste0(main_dir,'/cov/',geneid,'.txt'),sep='\t',quote = F,row.names = F)
          #output weights
          write.table(weights_df,paste0(main_dir,'/weights/',geneid,'.txt'),sep='\t',quote = F,row.names = F)
        }
        
      }else{
        print('INFO failed to build the model')
      }
    }else{
      print('INFO only one SNP available for the gene')
    }
  }else{
    print('INFO no genotype data for this gene')
  }
}


#------------------------------------------------

weights_cov<-function(main_dir,plink_file_name,expression_file_name,output_file_name,annotation_file_name){
  
  #mkdir
  if(!dir.exists(paste0(main_dir,'/output/'))){
    dir.create(paste0(main_dir,'/output/'))
  }
  
  #load annotation
  anno<-read.table(paste0(annotation_file_name),header = T,stringsAsFactors = F)
  
  #gene list
  gene_list<-sub('....$','',dir(paste0(main_dir,'/cov/')))
  
  #load weights
  print('INFO loading weights...')
  print(paste0('INFO totally ',length(gene_list),' imputable genes'))
  db_weights<-list()
  n_rows=0
  for (i in 1:length(gene_list)){
    if(i %in% c(seq(1,1000)*100,length(gene_list))){
      print(paste0('     ',i,' loaded'))
    }
    geneid<-gene_list[i]
    db_weights[[i]]<-read.table(paste0(main_dir,'/weights/',geneid,'.txt'),header = T,stringsAsFactors = F)
    n_rows=n_rows+nrow(db_weights[[i]])
    
    #cleaning tmp files
    if(opt$clean_tmp){
      file.remove(paste0(main_dir,'/weights/',geneid,'.txt'))
    }
    
  }
  
  #db.weight
  print('INFO generating db file...')
  weights<-as.data.frame(matrix(data=NA,ncol=4,nrow=n_rows))
  colnames(weights)<-c('rsid','gene','weight','eff_allele')
  
  #db.sample_info
  n_sample<-read.table(paste0(expression_file_name),header = T,nrows=1,stringsAsFactors = F)
  sample_info<-data.frame('n.samples'=ncol(n_sample)-1)
  
  #db.extra
  extra<-as.data.frame(matrix(data=NA,ncol=6,nrow=length(gene_list)))
  colnames(extra)<-c('gene','genename','pred.perf.R2','n.snps.in.model','pred.perf.pval','pred.perf.qval') 
  
  n_rows=0
  for (i in 1:length(gene_list)){
    #extra.gene
    extra[i,1]<-gene_list[i]
    #extra.genename
    extra[i,2]<-anno[which(anno$geneid==gene_list[i]),'genename']
    #extra.r2
    extra[i,3]<-db_weights[[i]][1,3]
    #extra.n.snps
    extra[i,4]<-nrow(db_weights[[i]])
    #extra.n.snps
    extra[i,5:6]<-db_weights[[i]][1,4] 
    
    weights[(n_rows+1):(n_rows+nrow(db_weights[[i]])),]<-as.matrix(data.frame('rsid'=db_weights[[i]]$rsid,'gene'=gene_list[i],'weight'=db_weights[[i]]$weights,'effect_allele'=db_weights[[i]]$effect_allele))
    n_rows=n_rows+nrow(db_weights[[i]])
    
  }
  
  extra$pred.perf.qval=p.adjust(extra$pred.perf.pval,method = 'BH')
  
  #load bim to get ref allele
  bim<-read.table(paste0(plink_file_name,'.bim'),header = F,stringsAsFactors = F)
  bim<-bim[bim[,2]!='.',];bim<-bim[!duplicated(bim[,2]),]
  bim<-bim[which(bim[,2] %in% weights[,1]),-1]
  
  weights<-merge(weights,bim,by=1)
  weights$ref_allele<-ifelse(weights$eff_allele==weights$V5,weights$V6,weights$V5)
  weights<-weights[,c(1,2,3,9,4)]
  
  #output db file
  print('INFO write db file...')
  
  if(file.exists(paste0(main_dir,"/output/",output_file_name,".db"))){
    file.remove(paste0(main_dir,"/output/",output_file_name,".db"))
  }
  
  db<-dbConnect(RSQLite::SQLite(), paste0(main_dir,"/output/",output_file_name,".db"))
  dbWriteTable(db, "weights", weights)
  dbWriteTable(db, "extra", extra)
  dbWriteTable(db, "sample_info", sample_info)
  dbDisconnect(db)
  
  #load cov
  print('INFO loading cov files...')
  cov<-list()
  n_rows=0
  for (i in 1:length(gene_list)){
    if(i %in% c(seq(1,1000)*100,length(gene_list))){
      print(paste0('     ',i,' loaded'))
    }
    geneid<-gene_list[i]
    cov[[i]]<-read.table(paste0(main_dir,'/cov/',geneid,'.txt'),header = T,stringsAsFactors = F)
    n_rows=n_rows+nrow(cov[[i]])
    
    #cleaning tmp files
    if(opt$clean_tmp){
      file.remove(paste0(main_dir,'/cov/',geneid,'.txt'))
    }
    
  }
  cov_df<-as.data.frame(matrix(data=NA,ncol=4,nrow=n_rows)) 
  colnames(cov_df)<-c('GENE','RSID1','RSID2','VALUE')
  #combine cov
  print('INFO combining cov file...')
  n_rows=0
  for (i in 1:length(gene_list)){
    cov_df[(n_rows+1):(n_rows+nrow(cov[[i]])),]<-cov[[i]]
    n_rows=n_rows+nrow(cov[[i]])
  }
  
  #output cov
  print('INFO writing cov file...')
  write.table(cov_df,paste0(main_dir,'/output/',output_file_name,'.cov'),sep='\t',row.names = F,quote = F)
  
  print('INFO finished')
}


#asso test func
asso<-function(gene_id){
  
  df = weights[weights$gene == gene_id, ]
  df$weight = as.numeric(df$weight)
  df$id = seq(1,nrow(df))
  
  cov_tmp = cov_all[cov_all$GENE == gene_id, ]
  cov_tmp = merge(cov_tmp, df[,c('rsid','id')], by.x = 'RSID1', by.y = 'rsid')
  cov_tmp = merge(cov_tmp, df[,c('rsid','id')], by.x = 'RSID2', by.y = 'rsid')
  
  cov_matrix = matrix(data = NA, nrow = nrow(df), ncol = nrow(df))
  
  for(i in 1:nrow(cov_tmp)){
    cov_matrix[cov_tmp[i,'id.x'], cov_tmp[i,'id.y']] = cov_tmp[i, 'VALUE']
    cov_matrix[cov_tmp[i,'id.y'], cov_tmp[i,'id.x']] = cov_tmp[i, 'VALUE']
  }
  
  #var(GReX)
  var_GReX = as.numeric(t(df$weight) %*% cov_matrix %*% df$weight)
  
  #beta
  effect_size = sum(df$weight * df[,gwas_beta_col] * diag(cov_matrix) / var_GReX)
  
  #z-score
  zscore = sum(df$weight * sqrt(diag(cov_matrix)) * df[,gwas_beta_col] / sqrt(var_GReX) / df[,gwas_se_col], na.rm = T)
  ans = list()
  ans[['effect_size']] = effect_size
  ans[['zscore']] = zscore
  ans[['n_snps_used']] = nrow(df)
  return(ans)
  
}


#model training
if(opt$model_training){
  print('INFO model training')
  print(paste0('INFO subjob id: ',opt$subjob_id))
  
  model_training(main_dir=opt$main_dir,
                 cis_window_size=opt$cis_window_size,
                 plink_file_name=opt$plink_file_name,
                 expression_file_name=opt$expression_file_name,
                 annotation_file_name=opt$annotation_file_name,
                 i=opt$subjob_id,
                 n_genes_for_each_subjob=opt$n_genes_for_each_subjob
  )
}


#generate db cov
if(opt$generate_db_and_cov){
  print('INFO generating db and cov files')
  
  weights_cov(main_dir=opt$main_dir,
              plink_file_name=opt$plink_file_name,
              expression_file_name=opt$expression_file_name,
              annotation_file_name=opt$annotation_file_name,
              output_file_name=opt$output_file_name
  )
}




#association test
if(opt$asso_test){
  db_path = opt$db_path
  cov_path = opt$cov_path
  gwas_path = opt$gwas_path
  gwas_variant_col = opt$gwas_variant_col
  gwas_beta_col = opt$gwas_beta_col
  gwas_or_col = opt$gwas_or_col
  gwas_se_col = opt$gwas_se_col
  gwas_p_col = opt$gwas_p_col
  gwas_zscore_col = opt$gwas_zscore_col
  gwas_eff_allele_col = opt$gwas_eff_allele_col
  gwas_ref_allele_col = opt$gwas_ref_allele_col
  gene_list_path = opt$gene_list_path_for_asso_test
  asso_out_path = opt$asso_out_path
  parallel = opt$parallel
  
  #load data
  cat('INFO loading weight file ...\n')
  con <- dbConnect(RSQLite::SQLite(), dbname=db_path) #establish connections
  weights = dbReadTable(con,"weights")
  colnames(weights) = c('rsid', 'gene', 'weight', 'model_ref_allele', 'model_eff_allele')
  extra_info = dbReadTable(con,"extra")
  dbDisconnect(con) #disconnect
  
  cat('INFO loading covariance ...\n')
  cov_all = data.frame(fread(cov_path))
  
  cat('INFO loading GWAS results ...\n')
  gwas <- data.frame(fread(gwas_path))

  #convert beta, or, se, pval
  gwas <- gwas[!duplicated(gwas[,gwas_variant_col]),]
  if(is.na(gwas_beta_col)){
    gwas$beta = log(gwas[,gwas_or_col],2.718)
    gwas_beta_col = 'beta'
  }
  
  if(is.na(gwas_se_col)){
    if(is.na(gwas_zscore_col)){
      gwas$z_score = qnorm(1-gwas[,gwas_p_col]/2)
      gwas$se = abs(gwas[,gwas_beta_col] / gwas[,'z_score'])
      gwas_se_col = 'se'
      
    }else{
      gwas$se = gwas[,gwas_beta_col] / gwas[,gwas_zscore_col]
      gwas_se_col = 'se'
    }
  }
  
  #upper allele
  gwas[,gwas_eff_allele_col] = toupper(gwas[,gwas_eff_allele_col])
  gwas[,gwas_ref_allele_col] = toupper(gwas[,gwas_ref_allele_col])
  
  #get the list of snps both available from gwas and prediction model
  weights <- merge(weights, gwas, by.x = 'rsid', by.y = gwas_variant_col)
  
  #strand check
  cat('INFO strand check and allele flipping ...\n')
  weights$allele_idx1 = weights[,'model_eff_allele'] == weights[,gwas_eff_allele_col]
  weights$allele_idx2 = weights[,'model_eff_allele'] == weights[,gwas_ref_allele_col]
  
  same_strand_pos = which(weights$allele_idx1 | weights$allele_idx2)
  num_of_variants_have_same_strand = length(same_strand_pos)
  
  if (num_of_variants_have_same_strand / nrow(weights) <0.9){
    stop('More than 10% of the variants in GWAS summary statistics have different strand, please check')
  }
  
  weights = weights[same_strand_pos,]
  
  #flip alleles
  weights[,gwas_beta_col] = ifelse(weights$allele_idx1, weights[,gwas_beta_col], weights[,gwas_beta_col]*-1)
  
  
  gene_list <- unique(weights$gene)
  
  if(!is.na(gene_list_path)){
    specified_gene_list = read.table(gene_list_path, header = T,stringsAsFactors = F)
    gene_list = intersect(gene_list, specified_gene_list[,1])
  }
  
  
  cat(paste0('INFO Association test for ',length(gene_list),' imputable genes ...\n'))
  
  
  if(parallel){
    registerDoParallel(cores = max(detectCores()-1, 1))
    
    trials <- length(gene_list)
    
    #system.time({
      result <- foreach(i = 1:trials, .combine=rbind) %dopar% {
        gene_id = gene_list[i]
        asso_ans = asso(gene_id)
        cat(paste0('INFO ',gene_id,' completed \n'))
        unlist(asso_ans)
      }
    #})
    
    stopImplicitCluster()
    
  }else{
    #The non-parallele way
    result <- as.data.frame(t(sapply(gene_list, function(x) asso(x))))
    result$effect_size = as.numeric(result$effect_size)
    result$zscore = as.numeric(result$zscore)
    result$n_snps_used = as.numeric(result$n_snps_used)
  }
  
  #output
  result = as.data.frame(result, stringsAsFactors = F)
  result$gene = gene_list
  result$pvalue = (1-pnorm(abs(result$zscore)))*2  #two sided
  result = result[which(result$effect_size * result$zscore>0),]
  result = merge(result, extra_info, by='gene')
  
  result = result[,c('gene','genename','zscore','effect_size','pvalue','pred.perf.R2','pred.perf.pval','n_snps_used','n.snps.in.model')]
  colnames(result) = c('gene','genename','zscore','effect_size','pvalue','pred_perf_R2','pred_perf_pval','n_snps_used','n_snps_in_model')
  
  result = result[order(abs(result$zscore),decreasing=T),]
  write.table(result, asso_out_path, quote = F, row.names = F, sep = '\t')
  
  cat('INFO done\n')
  
  
}












