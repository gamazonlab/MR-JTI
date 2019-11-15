#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

double vec_inner_prod(double* XX, double* theta, int t, int j, int M)
{
	double ans = 0.0;
	//theta[t*M+j] = 0.0;
	for(int k=0; k<M; k++){
		ans += XX[j*M+k]*theta[t*M+k];
	}
	return(ans);
}

double beta_j_norm(double* vec, int P)
{
	double ans = 0;
	for(int t=0; t<P; t++){
		ans += vec[t]*vec[t];
	}
	return(sqrt(ans));
}

void inner_iter(SEXP R_XX, double* XY, double* theta, int M, int P, double* beta_j_lasso, double* lambda1, double lambda2, double* Xnorm)
{
	double* tmp_XX;
	double tmp_XYj = 0.0;
	double bnorm = 0.0;
	for(int j=0; j<M; j++){
		for(int t=0; t<P; t++){
			theta[t*M+j] = 0.0;
			tmp_XX = REAL(VECTOR_ELT(R_XX, t));
			tmp_XYj = XY[t*M+j] - vec_inner_prod(tmp_XX, theta, t, j, M);
			beta_j_lasso[t] = fmax(fabs(tmp_XYj) - lambda1[t], 0)*copysign(1.0, tmp_XYj);
		}
		bnorm = beta_j_norm(beta_j_lasso, P);
		if(bnorm > 0){
			for(int t=0; t<P; t++){
				theta[t*M+j] = fmax(1 - lambda2/bnorm, 0)*beta_j_lasso[t]/Xnorm[j];
			}
		}
	}
}

SEXP wrapper(SEXP R_XX, SEXP R_XY, SEXP R_theta, SEXP R_M, SEXP R_P, SEXP R_beta_j_lasso, SEXP R_lambda1, SEXP R_lambda2, SEXP R_Xnorm)
{
	SEXP answer;
	double* ans;
	double* XY;
	double* theta;
	double* beta_j_lasso;
	double* Xnorm;
	int M, P;
	double* lambda1;
	double lambda2;
	XY = REAL(R_XY);
	theta = REAL(R_theta);
	beta_j_lasso = REAL(R_beta_j_lasso);
	Xnorm = REAL(R_Xnorm);
	M = INTEGER(R_M)[0];
	P = INTEGER(R_P)[0];
	lambda1 = REAL(R_lambda1);
	lambda2 = REAL(R_lambda2)[0];
	inner_iter(R_XX, XY, theta, M, P, beta_j_lasso, lambda1, lambda2, Xnorm);
	
	PROTECT(answer = allocVector(REALSXP,1));
	REAL(answer)[0] = 1;
  	UNPROTECT(1);
  	return(answer);
}
