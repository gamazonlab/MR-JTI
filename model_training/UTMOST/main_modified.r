
#cross tissue prediction UTMOST *modified
#gtex v8

args=as.numeric(commandArgs(TRUE))
print(args)
options(stringsAsFactors=F)
library(glmnet)
library(foreach)

#load glasso function
source('/home/zhoud2/script/v8/utmost/glasso.r')

#gene list
gene_list_geno<-dir('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_1m/')
gene_list_exp<-sub('....$','',dir('/data/coxvgi/zhoud2/projects/cross_tissue/exp_gene/'))
gene_list<-intersect(gene_list_geno,gene_list_exp)

#10 genes per subjob
i_start=min((args-1)*10+1,length(gene_list))
i_end=min(args*10,length(gene_list))

#gene_id='ENSG00000134243' #SORT1
#gene_id='ENSG00000198670' #LPA

for (gene_i in i_start:i_end){
  gene_id=gene_list[gene_i]
  print(gene_id)
  
  #mk output dir for each gene
  dir.create(paste0('/data/coxvgi/zhoud2/projects/cross_tissue/utmost/exp/',gene_id))
  
  #format transfer
  #expression
  exp<-readRDS(paste0('/data/coxvgi/zhoud2/projects/cross_tissue/exp_gene/',gene_id,'.rds'))
  tissue_list<-unique(exp$tissue)
  for (tissue in tissue_list){
    exp_tmp<-exp[exp$tissue==tissue,-3]
    saveRDS(exp_tmp,paste0('/data/coxvgi/zhoud2/projects/cross_tissue/utmost/exp/',gene_id,'/',tissue,'.rds'))
  }
  
  #genotype dosage
  geno<-read.table(paste0('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_1m/',gene_id),header = T,stringsAsFactors = F)
  geno[,1]<-paste0('GTEX.',geno[,1])
  saveRDS(geno,paste0('/data/coxvgi/zhoud2/projects/cross_tissue/utmost/geno/',gene_id,'.rds'))
  
  #genotype info
  geno_info<-read.table(paste0('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_1m_info/',gene_id),header = T,stringsAsFactors = F)
  chr_pos<-geno_info[,c(1,2)]
  geno_info<-geno_info[,-2]
  colnames(geno_info)<-c('SNP','REF.0.','ALT.1.')
  write.table(geno_info,paste0('/data/coxvgi/zhoud2/projects/cross_tissue/utmost/info/',gene_id),quote = F,row.names = F,sep='\t')
  
  ### data import ###
  main_dir='/data/coxvgi/zhoud2/projects/cross_tissue/utmost/'
  
  dose_path = paste0(main_dir,'geno/',gene_id,'.rds') #dosage
  info_path = paste0(main_dir,'info/',gene_id) #geno info
  Yt_path=paste0(main_dir,'exp/',gene_id,'/')  #expression
  Yt = dir(paste0(main_dir,'exp/',gene_id,'/'))  #tissue list
  ntune = 5  #number of grids for each tuning parameter
  
  P = length(Yt)
  fold = 5
  if(P){
    ## expr files ##
    Y = list()
    for(t in 1:P){
      Y[[t]] = readRDS(paste0(Yt_path, '/', Yt[t]))
    }
    ssize = unlist(lapply(Y, nrow))  #N of samples for each tissue
    T_num = length(Yt)  #N of tissues
    
    ## genotype files ##
    dose = readRDS(dose_path)
    dose$id<-dose[,1]
    dose<-dose[,c(1,ncol(dose),2:(ncol(dose)-1))]
    
    for(j in 3:ncol(dose)){ #center dosage to 0, not sure why
      dose[,j] = dose[,j] - mean(dose[,j])
    }
    N = nrow(dose) #overall genotyped sample size i.e. 838
    
    ## covariance matrix ##
    tmp = as.matrix(dose[,-(1:2)])
    XX = t(tmp)%*%as.matrix(tmp)/N #var cor matrix
    Xnorm = diag(XX)  #var of dosage
    remove(tmp); remove(XX)
    sub_id = dose[,1]  #subject id 838
    M = ncol(dose) - 2 #no of snps
    
    sub_id_map = list()  #subject id map for each tissue
    for(t in 1:T_num){
      tmp = rep(0, nrow(Y[[t]]))
      for(j in 1:length(tmp)){
        tmp[j] = which(sub_id==Y[[t]][j,1])
      }
      sub_id_map[[t]] = tmp
    }
    
    cv_config = cv_helper(N, fold) #cv setting
    cv_perm = cv_config$perm
    cv_idx = cv_config$idx
    
    single_res_test = list()
    single_lam = matrix(0,fold,P)
    single_theta_est = list()
    
    multi_res_test = list()
    multi_lam = matrix(0,fold,2)
    multi_theta_est = list()
    
    multi_res_test2 = list()
    multi_lam2 = array(0, dim=c(fold, P, 2))
    multi_theta_est2 = list()
    
    res_tune = list()
    rec_lamv = matrix(0, fold, ntune)
    
    avg_tune_res<-list() #average tuning error
    single_initial_est<-list()
    
    
    #----------------------
    
    for(f in 1:fold){
      print(f)
      #bgt = Sys.time()
      test_index = cv_perm[cv_idx[f,1]:cv_idx[f,2]]
      test_id = sub_id[test_index] #sample id in the testing set
      tuning_index = cv_perm[cv_idx[f%%fold+1,1]:cv_idx[f%%fold+1,2]]
      tuning_id = sub_id[tuning_index] #sample in the 'early stop' set
      
      #dfs for each set
      X_test = list()
      Y_test = list()
      X_tune = list()
      Y_tune = list()
      X_train = list()
      Y_train = list()
      X_all = list()
      Y_all = list()
      
      for(t in 1:T_num){
        X_train_tmp = sub_id_map[[t]][!(sub_id_map[[t]]%in%c(test_index,tuning_index))]
        Y_train_tmp = !(sub_id_map[[t]]%in%c(test_index,tuning_index))
        X_tuning_tmp = sub_id_map[[t]][(sub_id_map[[t]]%in%tuning_index)]
        Y_tuning_tmp = (sub_id_map[[t]]%in%tuning_index)
        X_test_tmp = sub_id_map[[t]][(sub_id_map[[t]]%in%test_index)]
        Y_test_tmp = (sub_id_map[[t]]%in%test_index)
        X_train[[t]] = apply(as.matrix(dose[X_train_tmp,-c(1,2)]),2,as.numeric)
        Y_train[[t]] = Y[[t]][Y_train_tmp, 2]
        X_tune[[t]] = apply(as.matrix(dose[X_tuning_tmp,-c(1,2)]),2,as.numeric)
        Y_tune[[t]] = Y[[t]][Y_tuning_tmp, 2]
        X_test[[t]] = apply(as.matrix(dose[X_test_tmp,-c(1,2)]),2,as.numeric)
        Y_test[[t]] = Y[[t]][Y_test_tmp, 2]
        X_all_tmp = sub_id_map[[t]]
        X_all[[t]] = apply(as.matrix(dose[X_all_tmp,-c(1,2)]),2,as.numeric)
        Y_all[[t]] = Y[[t]][,2]
        
        #rbind training and test set
        X_train[[t]]=rbind(X_train[[t]],X_test[[t]])
        Y_train[[t]]=c(Y_train[[t]],Y_test[[t]])
        
      }
      
      
      ## model training using all of the samples to get lambda range and initial est ##
      if (f==1){
        single_initial_est_all = matrix(0, ncol(X_train[[1]]), T_num) #N of snp by N tissue matrix
        single_summary_all = list()
        for(t in 1:T_num){
          set.seed(as.numeric(sub('^....','',gene_id)))
          tt = cv.glmnet(X_all[[t]], Y_all[[t]], alpha = 0.5, nfolds = 5,nlambda=50,pmax=200) #cv in 100% of the data, do not find the lambda here
          single_summary_all[[t]] = tt
          single_initial_est_all[,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)] #single tissue EN, beta and lambda with best cv error
        }
      }
      
      ## for each fold
      single_initial_est[[f]] = matrix(0, ncol(X_train[[1]]), T_num) #N of snp by N tissue matrix
      #single_summary = list()
      for(t in 1:T_num){
        set.seed(f)
        #print(length(Y_tune[[t]]))
        tt = cv.glmnet(X_train[[t]], Y_train[[t]], alpha = 0.5, nfolds = 5,nlambda=50,pmax=200) #cv in 100% of the data, do not find the lambda here
        #single_summary[[t]] = tt
        single_initial_est[[f]][,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)] #single tissue EN, beta and lambda with best cv error
      }
      
      
      
      ## use elastic net ests row norm as weights ## 
      lam_range = minmax_lambda(single_summary_all) #lambda range
      sig_norm = apply(single_initial_est[[f]], 1, function(x){sqrt(sum(x^2))}) #norm for each snp
      sig_norm[sig_norm==0] = rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2  #replace =0 with min?
      sig_norm = sig_norm/sum(sig_norm) #??? length=N of snps
      weights2 = 1/sig_norm; weights2 = weights2/sum(weights2);
      
      tis_norm = apply(single_initial_est[[f]], 2, function(x){sum(abs(x))})
      tis_norm[tis_norm==0] = rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
      tis_norm = tis_norm/sum(tis_norm)
      weights1 = 1/tis_norm; weights1 = weights1/sum(weights1);
      lam_V = seq(lam_range[1], lam_range[2], length.out = ntune)
      
      initial_numeric = as.numeric(single_initial_est[[f]])
      
      ## preparation
      XY = grad_prep(X_train, Y_train)
      XX_train = lapply(X_train, function(x){t(x)%*%x/nrow(x)})
      spsz = unlist(lapply(X_train,nrow)) #N of samples for each training set
      res_tune[[f]] = array(NA, dim=c(ntune, ntune, P)) #5*5 matrix *49 tissue, to find the best comb of lam1 and lam2
      
      rec_lamv[f,] = lam_V  #should be the same across 5 folds
      
      for(lam1 in 1:ntune){
        for(lam2 in 1:ntune){
          print(paste0('lam1= ',lam1,' lam2= ',lam2))
          single_est = matrix(initial_numeric, M, P) #snps by tissue matrix. initial iteration with single tissue estimates, save compute time
          ans = glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[lam1]/spsz, lambda2=lam_V[lam2], theta=single_est)
          
          
          #apply to tuning set
          if (lam1==1 & lam2==1){
            multi_res_test[[f]]<-list()
          }
          multi_res_test[[f]][[paste0(lam1,'_',lam2)]] = multi_mse(ans$est, X_tune, Y_tune)
        
          if(sum(ans$est!=0)>0){
            res_tune[[f]][lam1,lam2, ] = ans$tune_err #tuning err for each lam pair
            remove(single_est); remove(ans);
          }else{
            res_tune[[f]][lam1,lam2, ] = ans$tune_err
            remove(single_est); remove(ans);
            #break
          }			
        }
      }
      
      avg_tune_res[[f]] = apply(res_tune[[f]], c(1,2), mean) #average tuning error across the tissues for each lambda pair
    }
    #print(rec_lamv)  #should be the same across 5 folds
    
    #-----------------------------
    #find the best comb of lambda1 and lambda2 across 5 folds
    avg_tune_cross_fold<-matrix(rep(1,ntune*ntune),nrow = ntune,ncol=ntune)
    for (i_tune in 1:ntune){
      for (j_tune in 1:ntune){
        avg_tune_cross_fold[i_tune,j_tune]=mean(c(avg_tune_res[[1]][i_tune,j_tune],avg_tune_res[[2]][i_tune,j_tune],avg_tune_res[[3]][i_tune,j_tune],avg_tune_res[[4]][i_tune,j_tune],avg_tune_res[[5]][i_tune,j_tune]),na.rm = T)
        avg_tune_cross_fold[i_tune,j_tune]<-ifelse(is.na(avg_tune_cross_fold[i_tune,j_tune]),1,avg_tune_cross_fold[i_tune,j_tune])
      }
    }
    
    #find the best comb for each fold
    best.lam = which(avg_tune_cross_fold == min(avg_tune_cross_fold), arr.ind = TRUE)[1,] 
    
    #predicted expression in tuning set under best lambda comb
    multi_res<-list()
    for(f in 1:ntune){
      multi_res[[f]]<-multi_res_test[[f]][[paste0(as.numeric(best.lam)[1],'_',as.numeric(best.lam)[2])]]
    }
    
    #combine 5 fold results
    cv_df<-list()
    cv_r<-cv_p<-c()
    for (tissue_i in 1:T_num){
      cv_df[[tissue_i]]<-multi_res[[1]][tissue_i][[1]]
      
      for(fold_i in 2:fold){
        cv_df[[tissue_i]]<-rbind(cv_df[[tissue_i]],multi_res[[fold_i]][tissue_i][[1]])
      }
      print(nrow(cv_df[[tissue_i]]))
      fit<-cor.test(cv_df[[tissue_i]][,1],cv_df[[tissue_i]][,2])
      cv_r[tissue_i]<-as.numeric(fit$estimate)
      cv_p[tissue_i]<-as.numeric(fit$p.value)
      
      cv_r[tissue_i]<-ifelse(is.na(cv_r[tissue_i]),0,cv_r[tissue_i])
      cv_p[tissue_i]<-ifelse(is.na(cv_p[tissue_i]),1,cv_p[tissue_i])
    }
    
    
    ## generate an estimate with whole data ##
    XY = grad_prep(X_all, Y_all)
    XX_all = lapply(X_all, function(x){t(x)%*%x/nrow(x)})
    ans = glasso_no_early_stopping(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=lam_V[best.lam[1]]/spsz, lambda2=lam_V[best.lam[2]], theta=single_initial_est_all)
    
    info = read.table(info_path, header=T, sep='\t')
    downstream_est = data.frame(info[,1:3], ans$est)
    colnames(downstream_est)[4:ncol(downstream_est)]<-Yt
    
    #-------------------------------------------------------------
    #calculate r and p (retrained)
    pred_exp_list<-list()
    for (tissue_i in 1:length(Yt)){
      tissue_name=sub('....$','',Yt[tissue_i])
      pred_exp_list[[tissue_i]]<-t(apply(X_all[[tissue_i]],function(x) x*downstream_est[,tissue_i+3],MARGIN = 1))
      cor_result<-cor.test(rowSums(pred_exp_list[[tissue_i]]),Y_all[[tissue_i]]) #retrained
      
      weight_df<-data.frame(gene=gene_id,rsid=downstream_est$SNP,chr_bp=NA,ref_allele=downstream_est$REF.0.,counted_allele=downstream_est$ALT.1.,weight=downstream_est[,tissue_i+3],r2=round((cv_r[tissue_i])^2,3),p=cv_p[tissue_i],lambda=paste0(round(lam_V[best.lam[1]],2),';',round(lam_V[best.lam[2]],2)),r2_rt=as.numeric(cor_result$estimate^2),p_rt=cor_result$p.value,iter=ans$iter)
      
      #rm 0
      weight_df<-weight_df[weight_df$weight!=0,]
      
      #anno with chr and pos
      weight_df$chr_bp<-sapply(weight_df$rsid,function(x) chr_pos[which(x==chr_pos[,1]),2])
      
      if(nrow(weight_df)>0 & cv_r[tissue_i]>0.1 & cv_p[tissue_i]<0.05){
        print(tissue_i)
        saveRDS(weight_df,paste0('/data/coxvgi/zhoud2/projects/v8/weights_cv/ut/',tissue_name,'/',gene_id,'.rds'))
      }
      
    }
    
    
  }
  
}



