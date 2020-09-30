#-----MR-JTI-----

cat('-----MR-JTI----- \n')

options(digits = 4)

library("optparse")
library("glmnet")
library('HDCI')


option_list = list(
  make_option("--df_path", action="store", default=NA, type='character',
              help="Path to dataframe of gwas and eQTL summary statistics [required]"),
  make_option("--n_folds", action="store", default=5, type='integer',
              help="Number of cross-validation folds"),
  make_option("--n_snps", action="store", default=20, type='integer',
              help="The min number of SNPs (obs) to run MR-JTI"),
  make_option("--n_bootstrap", action="store", default=500, type='integer', 
              help="The number of bootstrap times"),
  make_option("--n_genes", action="store", default=1, type='integer',
              help="number of genes tested (for multiple testing correction)"),
  make_option("--weighted", action="store_true", default=FALSE,
              help="weighted by the inverse variance of GWAS effect size [default: %default]"),
  make_option("--out_path", action="store", default=NA, type='character',
              help="Path to output files [required]")
)

opt = parse_args(OptionParser(option_list=option_list))


#---Threshold Residual Bootstrap LASSO---
TRB_LASSO<-function (x, y, B = 500, alpha = 0.05, nfolds = 5, foldid, cv.OLS = FALSE, tau = 0, parallel = FALSE, standardize = TRUE, intercept = TRUE, parallel.boot = FALSE, ncores.boot = 1, weights, penalty.factor) {
  #   
  x <- as.matrix(x)
  y <- as.numeric(y)
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  selectset <- rep(0, p)
  Beta <- rep(0, p)
  
  #fit a global model
  globalfit <- glmnet(x, y, standardize = standardize, intercept = intercept, weights=weights, penalty.factor=penalty.factor)
  
  #get the lambda list
  lambda <- globalfit$lambda
  
  #cv glmnet
  cvfit <- weighted.escv.glmnet(x, y, lambda = lambda, nfolds = nfolds, tau = tau, cv.OLS = cv.OLS, standardize = standardize,intercept = intercept, weights=weights, penalty.factor=penalty.factor)
  
  #find the lambda with min cv.error
  lambda.opt <- cvfit$lambda.cv
  
  #beta including intercept
  fitlasso <- predict(globalfit, type = "coefficients",s = lambda.opt)
  #rm the intercept
  Beta <- fitlasso[-1] 
  
  #threshold beta
  Beta<-ifelse(abs(Beta)<1/ncol(x),0,Beta) 
  
  #predicted values
  fit_value <- predict(globalfit, newx = x, s = lambda.opt)
  
  #get residual
  residual <- y - fit_value
  residual_center <- residual - mean(residual)
  
  #beta matrix
  Beta.boot <- matrix(0, nrow = B, ncol = p)
  
  for (i in 1:B) {
    
    #resample bootstrap
    resam <- sample(1:n, n, replace = TRUE)
    ystar <- fit_value + residual_center[resam]
    
    #lasso
    boot.obj <- weighted.Lasso(x = x, y = ystar, lambda = lambda.opt, standardize = standardize, intercept = intercept, weights=weights, penalty.factor=penalty.factor)
    
    #get beta(alpha)
    Beta.boot[i,] <- boot.obj$beta
    
    #threshold beta(alpha)
    Beta.boot[i,]<-ifelse(abs(Beta.boot[i,])<1/ncol(x),0,Beta.boot[i,])
  }
  
  #mean of beta(alpha)
  beta_mean<-apply(Beta.boot,2,function(x) mean(x))
  
  #CI
  interval <- matrix(0, 2, p)
  bound.percentile <- apply(Beta.boot, 2, function(u) {
    quantile(u, prob = c(1 - alpha/2, alpha/2))
  })
  interval[1, ] <- 2 * Beta - bound.percentile[1, ]
  interval[2, ] <- 2 * Beta - bound.percentile[2, ]
  
  #output
  object <- list(lambda.opt = lambda.opt, beta = beta_mean, interval = interval)
  object
}


#estimation stability with cross-validation (based on 'HDCI')
weighted.escv.glmnet<-function (x, y, lambda = NULL, nfolds = 5, cv.OLS = FALSE, tau = 0, standardize = TRUE, intercept = TRUE, weights, penalty.factor) 
{
  
  if (!is.null(lambda) && length(lambda) < 2) {
    stop("Need more than one value of lambda for escv.glmnet")
  }
  n <- nrow(x)
  p <- ncol(x)
  y <- drop(y)
  glmnet.call <- match.call(expand.dots = TRUE)
  which <- match(c("nfolds", "foldid"), names(glmnet.call), 
                 F)
  if (any(which)) {
    glmnet.call <- glmnet.call[-which]
  }
  glmnet.call[[1]] <- as.name("glmnet")
  glmnet.object <- glmnet(x, y, lambda = lambda, standardize = standardize,intercept = intercept, weights = weights, penalty.factor=penalty.factor) #
  
  lambda <- glmnet.object$lambda
  glmnet.object$call <- glmnet.call
  
  foldid <- sample(rep(seq(nfolds), length = n))
  
  out <- list()
  
  for (k in 1:nfolds) {
    test <- foldid == k
    train <- foldid != k
    obj <- glmnet(x[train, , drop = FALSE], y[train], 
                  lambda = lambda, standardize = standardize, intercept = intercept, weights = weights[which(train)], penalty.factor=penalty.factor) #
    fitmat <- predict(obj, newx = x)
    predtest <- predict(obj, newx = x[test, , drop = FALSE])
    residmat <- apply((y[test] - predtest)^2, 2, mean)
    out[[k]] <- list(residmat = residmat, fitmat = fitmat)
  }
  
  residmat <- matrix(0, length(lambda), nfolds)
  residmates <- matrix(0, length(lambda), nfolds)
  fitmat <- array(0, dim = c(n, length(lambda), nfolds))
  for (k in 1:nfolds) {
    residmat[, k] <- out[[k]]$residmat
    fitmat[, , k] <- out[[k]]$fitmat
    out[[k]]$residmat <- NULL
    out[[k]]$fitmat <- NULL
  }
  meanfit <- apply(fitmat, c(1, 2), mean)
  meanfit2 <- apply(meanfit^2, 2, sum)
  for (k in 1:nfolds) {
    residmates[, k] <- apply((fitmat[, , k] - meanfit)^2, 
                             2, sum)/meanfit2
  }
  residmates[is.na(residmates)] <- Inf
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/nfolds)
  es <- apply(residmates, 1, mean)
  es.error <- sqrt(apply(residmates, 1, var)/nfolds)
  indcv <- which.min(cv)
  lambda.cv <- lambda[indcv]
  indcv0 <- indcv
  cv1se <- cv
  cv1se[cv <= (cv[indcv] + cv.error[indcv])] <- cv[indcv] + 
    cv.error[indcv]
  indcv1se <- which.min(cv1se)
  lambda.cv1se <- lambda[indcv1se]
  indescv <- which.min(es[1:indcv])
  lambda.escv <- lambda[indescv]
  if (cv.OLS) {
    out <- list()
    
    for (k in 1:nfolds) {
      test <- foldid == k
      train <- foldid != k
      obj <- glmnet(x[train, , drop = FALSE], y[train], 
                    lambda = lambda[1:indcv], standardize = standardize, 
                    intercept = intercept, weights = weights[which(train)],penalty.factor=penalty.factor) #
      fitmat <- predict(obj, newx = x)
      predtest <- predict(obj, newx = x[test, , drop = FALSE])
      selectset0 <- rep(0, p)
      for (i in 1:indcv) {
        selectset <- abs(obj$beta[, i]) > 0
        if (sum(selectset) > 0) {
          if (sum(abs(selectset - selectset0)) > 0) {
            mls.obj <- mls(x[train, selectset, drop = FALSE], 
                           y[train], tau = tau, standardize = standardize, 
                           intercept = intercept)
            fitmat[, i] <- mypredict(mls.obj, newx = x[, 
                                                       selectset, drop = FALSE])
            predtest[, i] <- mypredict(mls.obj, newx = x[test, 
                                                         selectset, drop = FALSE])
          }
          else {
            fitmat[, i] <- fitmat[, i - 1]
            predtest[, i] <- predtest[, i - 1]
          }
          selectset0 <- selectset
        }
      }
      residmat <- apply((y[test] - predtest)^2, 2, 
                        mean)
      out[[k]] <- list(residmat = residmat, fitmat = fitmat)
    }
    
    
    residmat <- matrix(0, indcv, nfolds)
    residmates <- matrix(0, indcv, nfolds)
    fitmat <- array(0, dim = c(n, indcv, nfolds))
    for (k in 1:nfolds) {
      residmat[, k] <- out[[k]]$residmat
      fitmat[, , k] <- out[[k]]$fitmat
      out[[k]]$residmat <- NULL
      out[[k]]$fitmat <- NULL
    }
    meanfit <- apply(fitmat, c(1, 2), mean)
    meanfit2 <- apply(meanfit^2, 2, sum)
    for (k in 1:nfolds) {
      residmates[, k] <- apply((fitmat[, , k] - meanfit)^2, 
                               2, sum)/meanfit2
    }
    residmates[is.na(residmates)] <- Inf
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/nfolds)
    es <- apply(residmates, 1, mean)
    es.error <- sqrt(apply(residmates, 1, var)/nfolds)
    indcv <- which.min(cv)
    lambda.cv <- lambda[indcv]
    cv1se <- cv
    cv1se[cv <= (cv[indcv] + cv.error[indcv])] <- cv[indcv] + 
      cv.error[indcv]
    indcv1se <- which.min(cv1se)
    lambda.cv1se <- lambda[indcv1se]
    indescv <- which.min(es[1:indcv0])
    lambda.escv <- lambda[indescv]
  }
  object <- list(lambda = lambda, glmnet.fit = glmnet.object, 
                 cv = cv, cv.error = cv.error, es = es, es.error = es.error, 
                 lambda.cv = lambda.cv, lambda.cv1se = lambda.cv1se, lambda.escv = lambda.escv)
  object
}



#Lasso (based on 'HDCI')
weighted.Lasso<-function (x, y, lambda = NULL, fix.lambda = TRUE, cv.method = "cv", nfolds = 5, foldid, cv.OLS = FALSE, tau = 0, standardize = TRUE, intercept = TRUE, weights, penalty.factor) 
{
  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (is.null(lambda)) {
    stop("Should give a value of lambda for fix.lambda=TRUE")
  }
  if (length(lambda) > 1) {
    stop("The length of lambda should be 1 if fix.lambda=TRUE")
  }
  
  globalfit <- glmnet(x, y, standardize = standardize, intercept = intercept,weights = weights, penalty.factor=penalty.factor) #
  fitlasso <- predict(globalfit, type = "coefficients", 
                      s = lambda)
  beta0 <- fitlasso[1]
  beta <- fitlasso[-1]
  
  
  if (intercept) {
    meanx <- apply(x, 2, mean)
    mu <- mean(y)
  }
  else {
    meanx <- rep(0, p)
    mu <- 0
  }
  object <- list()
  object$beta0 <- beta0
  object$beta <- beta
  object$lambda <- lambda
  object$meanx <- meanx
  object$mu <- mu
  object
}




#load df
cat('INFO loading dataframe... \nPlease make sure that the effect alleles of GWAS and eQTL are correctly aligned. \n')
df<-read.table(opt$df_path,header = T,stringsAsFactors = F)

#skip genes with less than 20(defult) obs (no way for cv)
if (nrow(df)<opt$n_snps){stop('Need >= 20 variant to run MR-JTI')}

#no of sig gwas loci (only consider these for pleiotropy control)
n_gwas<-length(which(df$gwas_p<0.05))

#order by gwas p value
df<-df[order(df$gwas_p),]

#weighted by the precision
weights=rep(1,nrow(df))
if(opt$weighted){
  weights=1/(df$gwas_se)^2
}

#penalty factor
penalty.factor<-c(0,0,rep(1,n_gwas))

#generate identity matrix
df<-df[,c('rsid','gwas_beta','eqtl_beta','ldscore')] 
if(n_gwas>0){
  df[,(ncol(df)+1):(ncol(df)+n_gwas)]<-diag(nrow(df))[,1:n_gwas]
}

#---assign x and y---
y=df[,'gwas_beta'];x=as.matrix(df[,3:ncol(df)])
#scale x and y
y=scale(y)
x=apply(x, 2, function(x) scale(x))

#---real data bootstrap lasso---
cat('INFO running residual bootstrap lasso... \n')
ans<-TRB_LASSO(x=x,y=y,intercept = T,standardize = T,nfolds = opt$n_folds,alpha=0.05/opt$n_genes,B=opt$n_bootstrap, weights=weights, penalty.factor=penalty.factor)

#---output---
output<-data.frame(variable=c('expression','ldsc',df[1:n_gwas,1]))
output$beta=ans$beta
output$beta_CI_lower=as.numeric(ans$interval[1,])
output$beta_CI_upper=as.numeric(ans$interval[2,])
output$CI_significance<-ifelse(output$beta_CI_lower*output$beta_CI_upper>0,'sig','nonsig')


#output
write.csv(output,opt$out_path,quote = F,row.names = F)

cat('INFO done \n-----\n')

