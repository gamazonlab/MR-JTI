
#cross tissue prediction UTMOST *original
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

#gene_id='ENSG00000134243'

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
    row.names(exp_tmp)<-seq(1,nrow(exp_tmp))
    saveRDS(exp_tmp,paste0('/data/coxvgi/zhoud2/projects/cross_tissue/utmost/exp/',gene_id,'/',tissue,'.rds'))
  }
  
  #genotype dosage
  geno<-read.table(paste0('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_1m/',gene_id),header = T,stringsAsFactors = F)
  geno[,1]<-paste0('GTEX.',geno[,1])
  row.names(geno)<-seq(1,nrow(geno))
  saveRDS(geno,paste0('/data/coxvgi/zhoud2/projects/cross_tissue/utmost/geno/',gene_id,'.rds'))
  
  #genotype info
  geno_info<-read.table(paste0('/data/coxvgi/zhoud2/data/gtex/geno/v8/gene/p_1m_info/',gene_id),header = T,stringsAsFactors = F)
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
    ssize = unlist(lapply(Y, nrow))
    T_num = length(Yt)
    
    ## genotype files ##
    dose = readRDS(dose_path)
    dose$id<-dose[,1]
    dose<-dose[,c(1,ncol(dose),2:(ncol(dose)-1))]
    
    for(j in 3:ncol(dose)){ ## if no 'dose' column
      dose[,j] = dose[,j] - mean(dose[,j])
    }
    N = nrow(dose)
    ## covariance matrix ##
    tmp = as.matrix(dose[,-(1:2)])
    XX = t(tmp)%*%as.matrix(tmp)/N
    Xnorm = diag(XX)
    remove(tmp); remove(XX)
    sub_id = dose[,1]
    M = ncol(dose) - 2
    sub_id_map = list()
    for(t in 1:T_num){
      tmp = rep(0, nrow(Y[[t]]))
      for(j in 1:length(tmp)){
        tmp[j] = which(sub_id==Y[[t]][j,1])
      }
      sub_id_map[[t]] = tmp
    }
    cv_config = cv_helper(N, fold)
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
    for(f in 1:fold){
      print(f)
      #bgt = Sys.time()
      test_index = cv_perm[cv_idx[f,1]:cv_idx[f,2]]
      test_id = sub_id[test_index]
      tuning_index = cv_perm[cv_idx[f%%fold+1,1]:cv_idx[f%%fold+1,2]]
      tuning_id = sub_id[tuning_index]
      
      X_test = list()
      Y_test = list()
      X_tune = list()
      Y_tune = list()
      X_train = list()
      Y_train = list()
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
      }
      
      ## model training ##	
      ## train elastic net and used average lambda as tuning parameters ##
      single_initial_est = matrix(0, ncol(X_train[[1]]), T_num)
      single_summary = list()
      for(t in 1:T_num){
        tt = cv.glmnet(X_train[[t]], Y_train[[t]], alpha = 0.5, nfolds = 5)
        single_summary[[t]] = tt
        single_initial_est[,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)]
      }
      ## performance of Elastic net on tuning and testing data with various tuning parameters
      els_output = elastic_net_mse(single_summary, X_tune, Y_tune, X_test, Y_test)
      single_res_test[[f]] = els_output$mse
      single_lam[f,] = els_output$lam
      single_theta_est[[f]] = els_output$est
      remove(els_output)
      
      ## use elastic net ests row norm as weights ##
      lam_range = minmax_lambda(single_summary)
      sig_norm = apply(single_initial_est, 1, function(x){sqrt(sum(x^2))})
      sig_norm[sig_norm==0] = rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2
      sig_norm = sig_norm/sum(sig_norm)
      weights2 = 1/sig_norm; weights2 = weights2/sum(weights2);
      
      tis_norm = apply(single_initial_est, 2, function(x){sum(abs(x))})
      tis_norm[tis_norm==0] = rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
      tis_norm = tis_norm/sum(tis_norm)
      weights1 = 1/tis_norm; weights1 = weights1/sum(weights1);
      lam_V = seq(lam_range[1], lam_range[2], length.out = ntune)
      #lam_V = seq(lam_range[1], lam_range[2], length.out = ntune)
      
      
      initial_numeric = as.numeric(single_initial_est)
      remove(single_summary); remove(single_initial_est);
      
      ## preparation
      XY = grad_prep(X_train, Y_train)
      XX_train = lapply(X_train, function(x){t(x)%*%x/nrow(x)})
      spsz = unlist(lapply(X_train,nrow))
      #res_tune = rep(0, ntune)
      res_tune[[f]] = array(-1, dim=c(ntune, ntune, P))
      #best.lam = 0
      rec_lamv[f,] = lam_V
      for(lam1 in 1:ntune){
        for(lam2 in 1:ntune){
          single_est = matrix(initial_numeric, M, P)
          ans = glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[lam1]/spsz, lambda2=lam_V[lam2], theta=single_est)
          if(sum(ans$est!=0)>0){
            res_tune[[f]][lam1,lam2, ] = ans$tune_err
            #cat("lambda1=",lam_V[lam1], "; lambda2=", lam_V[lam2], "; avg tune err=", ans$avg_tune_err, '\n')
            remove(single_est); remove(ans);
          }else{
            remove(single_est); remove(ans);
            break
          }			
        }
      }
      avg_tune_res = apply(res_tune[[f]], c(1,2), mean)
      best.lam = which(avg_tune_res == min(avg_tune_res[avg_tune_res>=0]), arr.ind = TRUE)[1,]
      single_est = matrix(initial_numeric, M, P)
      ans = glasso(X=X_train, Y=Y_train, X1=X_tune, Y1=Y_tune, XX=XX_train, XY=XY, Xnorm=Xnorm, lambda1=lam_V[best.lam[1]]/spsz, lambda2=lam_V[best.lam[2]], theta=single_est)
      multi_res_test[[f]] = multi_mse(ans$est, X_test, Y_test)
      multi_lam[f,] = lam_V[best.lam]
      multi_theta_est[[f]] = ans$est
      remove(single_est); remove(ans);

    }

    ## generate an estimate with whole data ##
    X_all = list()
    Y_all = list()
    for(t in 1:T_num){
      X_all_tmp = sub_id_map[[t]]
      X_all[[t]] = apply(as.matrix(dose[X_all_tmp,-c(1,2)]),2,as.numeric)
      Y_all[[t]] = Y[[t]][,2]
    }
    # initial values 
    single_initial_est = matrix(0, ncol(X_train[[1]]), T_num)
    for(t in 1:T_num){
      tt = cv.glmnet(X_all[[t]], Y_all[[t]], alpha = 0.5, nfolds = 5)
      single_initial_est[,t] = tt$glmnet.fit$beta[,which.min(tt$cvm)]
    }
    
    sig_norm = apply(single_initial_est, 1, function(x){sqrt(sum(x^2))})
    sig_norm[sig_norm==0] = rep(min(sig_norm[sig_norm>0]), sum(sig_norm==0))/2
    sig_norm = sig_norm/sum(sig_norm)
    weights2 = 1/sig_norm; weights2 = weights2/sum(weights2);
    
    tis_norm = apply(single_initial_est, 2, function(x){sum(abs(x))})
    tis_norm[tis_norm==0] = rep(min(tis_norm[tis_norm>0]), sum(tis_norm==0))/2
    tis_norm = tis_norm/sum(tis_norm)
    weights1 = 1/tis_norm; weights1 = weights1/sum(weights1);
    
    spsz = unlist(lapply(X_all,nrow))
    initial_numeric = as.numeric(single_initial_est)
    
    XY = grad_prep(X_all, Y_all)
    XX_all = lapply(X_all, function(x){t(x)%*%x/nrow(x)})
    tmp_res = rep(0, fold)
    for(f in 1:fold){
      ans = glasso_no_early_stopping(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=multi_lam[f,1]/spsz, lambda2=multi_lam[f,2], theta=matrix(initial_numeric,M,P))
      tmp_res[f] = ans$avg_train_err
    }
    final.lam = multi_lam[which.min(tmp_res),]
    ans = glasso_no_early_stopping(X=X_all, Y=Y_all, XX=XX_all, XY=XY, Xnorm=Xnorm, lambda1=final.lam[1]/spsz, lambda2=final.lam[2], theta=matrix(initial_numeric,M,P))
    info = read.table(info_path, header=T, sep='\t')
    downstream_est = data.frame(info[,1:3], ans$est)
    colnames(downstream_est)[4:ncol(downstream_est)]<-Yt

    #-------------------------------------------------------------
    #calculate r2 and p
    pred_exp_list<-list()
    for (tissue_i in 1:length(Yt)){
      tissue_name=sub('....$','',Yt[tissue_i])
      pred_exp_list[[tissue_i]]<-t(apply(X_all[[tissue_i]],function(x) x*downstream_est[,tissue_i+3],MARGIN = 1))
      cor_result<-cor.test(rowSums(pred_exp_list[[tissue_i]]),Y_all[[tissue_i]])
      weight_df<-data.frame(gene=gene_id,rsid=downstream_est$SNP,chr_bp=NA,ref_allele=downstream_est$REF.0.,counted_allele=downstream_est$ALT.1.,weight=downstream_est[,tissue_i+3],r2=as.numeric(cor_result$estimate^2),p=cor_result$p.value,lambda=paste0(round(final.lam[1],2),';',round(final.lam[2],2)),iter=ans$iter)
      weight_df<-weight_df[weight_df$weight!=0,]
      
      if(nrow(weight_df)>0 & cor_result$estimate>0.1 & cor_result$p.value<0.05){
        rownames(weight_df)<-seq(1,nrow(weight_df))
        saveRDS(weight_df,paste0('/data/coxvgi/zhoud2/projects/v8/weights_rt/ut/',tissue_name,'/',gene_id,'.rds'))
      }
      
    }
    
  }
  
}

















