
optim_path='/home/z***/script/v8/utmost/'

#args = commandArgs(trailingOnly=TRUE)
### optimization part ###
grad_prep <- function(X, Y){
  ## pre-calculate some metrics for gradient
  ## args
  ## X: a list of covariate matrices corresponding to each response
  ## Y: a list of response vectors
  ## value
  ## XY: a list of matrices X^TY for each response
  ll = length(Y)
  P = ncol(X[[1]])
  XY = matrix(0,P,ll)
  for(i in 1:ll){
    XY[,i] = t(X[[i]])%*%Y[[i]]/nrow(X[[i]])
  }
  XY
}

cv_helper <- function(N, fold){
  ## helper function for generating cross-validation sets
  ## args
  ## N: number of sample size
  ## fold: number of folds
  ## values
  ## perm: a permutation of 1 to N
  ## idx: matrix of fold by 2 with first col being starting index and second col being ending index
  valid_num = floor(N/fold)
  set.seed(123)
  perm = sample(1:N, size = N)
  idx1 = seq(1,N,valid_num)
  idx2 = c(idx1[-1]-1,N)
  list(perm=perm, idx=cbind(idx1,idx2))
}

minmax_lambda <- function(lst){
  ## get the minimum and maximum of lambda searched in cross-validation of an elastic net model
  ## args
  ## lst: an object returned by glmnet
  ## value
  ## min_lam: smallest lambda searched in glmnet cross-validation
  ## max_lam: largest lambda searched in glmnet cross-validation
  max_lam = max(unlist(lapply(lst, function(x){max(x$lambda)})))
  min_lam = min(unlist(lapply(lst, function(x){min(x$lambda)})))
  c(min_lam, max_lam)
}



elastic_net_mse <- function(lst, X_tune, Y_tune, X_test, Y_test){
  ## evaluate the performance of elastic net on each response
  ## args
  ## lst: a list of glmnet object (fitted elastic net model for each response)
  ## X_tune: a list of covariate matrices corresponding for each response (for tuning lambda)
  ## Y_tune: a list of response vectors (for tuning lambda)
  ## X_test: a list of covariate matrices corresponding for each response (for testing performance)
  ## Y_test: a list of response vectors (for testing performance)
  ## value
  ## lam: best performing lambda (on (X_tune,Y_tune)) for each response
  ## mse: list of matrices with each element being a matrix of predicted vs observed response
  ## est: estimated effect sizes for each response (B matrix)
  
  #lst<-single_summary
  P = length(lst)  #N of tissues
  M = ncol(X_tune[[1]]) #N of SNPs
  lam_V = rep(0, P)
  test_res = list()
  test_beta = matrix(0, M, P)
  for(t in 1:P){
    ncv = length(lst[[t]]$lambda)  #no of cv, i.e. no of tested lambdas
    tmp_mse = rep(0, ncv)
    for(k in 1:ncv){
      tmp_mse[k] = mean((Y_tune[[t]] - X_tune[[t]]%*%lst[[t]]$glmnet.fit$beta[,k])^2) #mean error in tuning set (apply beta to tuning set)
    }
    ss = which.min(tmp_mse)
    test_beta[,t] = lst[[t]]$glmnet.fit$beta[,ss]
    lam_V[t] = lst[[t]]$lambda[ss]
    predicted = X_test[[t]]%*%lst[[t]]$glmnet.fit$beta[,ss]
    test_res[[t]] = cbind(Y_test[[t]], predicted)
  }
  list(lam = lam_V, mse = test_res, est = test_beta)
}

multi_mse <- function(theta_est, X_test, Y_test){
  answer = list()
  P = ncol(theta_est)
  for(t in 1:P){
    predicted = X_test[[t]]%*%theta_est[,t]
    answer[[t]] = cbind(Y_test[[t]], predicted)
  }
  answer
}

avg_perm <- function(mse_lst){
  fd = length(mse_lst)
  P = length(mse_lst[[1]])
  rsq = mse = adj_mse = matrix(0, fd, P)
  for(f in 1:fd){
    for(t in 1:P){
      rsq[f,t] = (cor(mse_lst[[f]][[t]])[1,2])^2
      mse[f,t] = mean((mse_lst[[f]][[t]][,1]-mse_lst[[f]][[t]][,2])^2)
      adj_mse[f,t] = mse[f,t]/var(mse_lst[[f]][[t]][,1])
    }
  }
  cbind(apply(rsq, 2, mean), apply(mse, 2, mean), apply(adj_mse, 2, mean))
  
  #list(rsq = apply(rsq, 2, mean), mse = apply(mse, 2, mean), adj_mse = apply(adj_mse, 2, mean))
}

pred_test <- function(Y){
  if(sum(Y[,2]==0)==nrow(Y)|var(Y[,2])==0){
    return(2)
  }else{
    summary(lm(Y[,1]~Y[,2]))$coefficients[2,4]
  }
}

glasso <- function(X, Y, X1, Y1, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3){
  bgt = Sys.time()
  M = nrow(XY)
  P = length(X)
  NN = unlist(lapply(X, nrow))
  old_objV1 = rep(0,P)
  for(t in 1:P){
    old_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)  #average error, appling st weight on training set
  }
  old_objV2 = rep(0,P)
  for(t in 1:P){
    old_objV2[t] = 1/2*mean((Y1[[t]]-X1[[t]]%*%theta[,t])^2) #average error, appling st weight on tuning set
  }
  beta_j_lasso = rep(0, P)
  tmp_XYj = 0
  if(!is.loaded("wrapper")){
    dyn.load(paste0(optim_path,"optim.so"))
  }
  for(i in 1:maxiter){  #iterations to get the lowest error.
    bgt = Sys.time()
    res = .Call("wrapper", XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
    edt = Sys.time()
    #print(edt-bgt)
    new_objV1 = new_objV2 = rep(0,P)
    for(t in 1:P){
      new_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
    }
    #cat("Training error: ", new_objV1, '\n')
    for(t in 1:P){
      new_objV2[t] = 1/2*mean((Y1[[t]]-X1[[t]]%*%theta[,t])^2)
    }
    #cat("Testing error: ", new_objV2, '\n')
    if(mean(new_objV2) > mean(old_objV2)|mean(new_objV1) > mean(old_objV1)){
      break
    }else{
      old_objV2 = new_objV2
    }
    if(max(abs(new_objV1-old_objV1)) < eps){
      break
    }else{
      old_objV1 = new_objV1
    }
  }
  #print(paste0('iterations: ',i))
  #edt = Sys.time()
  #print(edt-bgt)
  list(est = theta, avg_tune_err = mean(new_objV2), tune_err=new_objV2)
}

glasso_no_early_stopping <- function(X, Y, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3){
  M = nrow(XY)
  P = length(X)
  NN = unlist(lapply(X, nrow))
  old_objV1 = rep(0,P)
  for(t in 1:P){
    old_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
  }
  cat("Training error: ", mean(old_objV1), '\n')
  beta_j_lasso = rep(0, P)
  tmp_XYj = 0
  if(!is.loaded("wrapper")){
    dyn.load(paste0(optim_path,"optim.so"))
  }
  for(i in 1:maxiter){
    res = .Call("wrapper", XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
    new_objV1 = rep(0,P)
    for(t in 1:P){
      new_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)
    }
    cat("Training error: ", mean(new_objV1), '\n')
    if(max(abs(new_objV1-old_objV1)) < eps|mean(new_objV1) > mean(old_objV1)){
      break
    }else{
      old_objV1 = new_objV1
    }
  }
  #print(paste0('iterations: ',i))
  list(est = theta, avg_train_err = mean(new_objV1), train_err = new_objV1,iter=i)
}

### command line input ###











glasso_mod <- function(X, Y, X1, Y1, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3){
  bgt = Sys.time()
  M = nrow(XY)
  P = length(X)
  NN = unlist(lapply(X, nrow))
  old_objV1 = rep(0,P)
  for(t in 1:P){
    old_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2)  #average error, appling st weight on training set
  }
  old_objV2 = rep(0,P)
  for(t in 1:P){
    old_objV2[t] = 1/2*mean((Y1[[t]]-X1[[t]]%*%theta[,t])^2) #average error, appling st weight on tuning set
  }
  beta_j_lasso = rep(0, P)
  tmp_XYj = 0
  if(!is.loaded("wrapper")){
    dyn.load(paste0(optim_path,"optim.so"))
  }
  
  beta_list<-list()
  beta_list[[1]]<-as.numeric(theta)
  
  for(i in 1:maxiter){  #iterations to get the lowest error.
    beta_matrix<-matrix(beta_list[[i]],M,P)
    
    bgt = Sys.time()
    res = .Call("wrapper", XX, XY, beta_matrix, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
    edt = Sys.time()
    #print(edt-bgt)
    new_objV1 = new_objV2 = rep(0,P)
    for(t in 1:P){
      new_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%beta_matrix[,t])^2)
    }
    #cat("Training error: ", new_objV1, '\n')
    for(t in 1:P){
      new_objV2[t] = 1/2*mean((Y1[[t]]-X1[[t]]%*%beta_matrix[,t])^2)
    }
    #cat("Testing error: ", new_objV2, '\n')
    if(mean(new_objV2) > mean(old_objV2)|mean(new_objV1) > mean(old_objV1)){
      beta_output<-beta_list[[i]]
      break
    }else{
      old_objV2 = new_objV2
      beta_list[[i+1]]<-as.numeric(beta_matrix)
    }
    if(max(abs(new_objV1-old_objV1)) < eps){
      beta_output<-beta_list[[i]]
      break
    }else{
      old_objV1 = new_objV1
      beta_list[[i+1]]<-as.numeric(beta_matrix)
    }
  }
  #print(paste0('iterations: ',i))
  #edt = Sys.time()
  #print(edt-bgt)
  list(est = matrix(beta_output,M,P), avg_tune_err = mean(new_objV2), tune_err=new_objV2)
}


glasso_no_early_stopping_mod <- function(X, Y, XX, XY, Xnorm, lambda1, lambda2, theta, stepsize = 1e-4, maxiter = 50, eps = 1e-3){
  M = nrow(XY)
  P = length(X) #n of tissues i.e. 49
  NN = unlist(lapply(X, nrow))
  old_objV1 = rep(0,P)
  for(t in 1:P){
    old_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%theta[,t])^2) #error of single tissue
  }
  cat("Training error: ", mean(old_objV1), '\n')
  beta_j_lasso = rep(0, P) #???
  tmp_XYj = 0
  if(!is.loaded("wrapper")){
    dyn.load(paste0(optim_path,"optim.so"))
  }
  
  beta_list<-list()
  beta_list[[1]]<-as.numeric(theta)
  
  for(i in 1:maxiter){
    
    beta_matrix<-matrix(beta_list[[i]],M,P)
    res = .Call("wrapper", XX, XY, beta_matrix, M, P, beta_j_lasso, lambda1, lambda2, Xnorm)
    
    new_objV1 = rep(0,P)
    for(t in 1:P){
      new_objV1[t] = 1/2*mean((Y[[t]]-X[[t]]%*%beta_matrix[,t])^2) #error after optimization
    }
    cat("Training error: ", mean(new_objV1), '\n')
    if(max(abs(new_objV1-old_objV1)) < eps|mean(new_objV1) > mean(old_objV1)){
      beta_output<-beta_list[[i]]
      break
    }else{
      old_objV1 = new_objV1
      beta_list[[i+1]]<-as.numeric(beta_matrix)
    }
  }
  #print(paste0('iterations: ',i))
  list(est = matrix(beta_output,M,P), avg_train_err = mean(new_objV1), train_err = new_objV1,iter=i)
}

### command line input ###











