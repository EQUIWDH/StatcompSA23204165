// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
//' @import Rcpp
//' @import RcppEigen
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;
//' @title Calculate group version Cp criterion
//' @description Use Rcpp for calculating Cp criterion
//' @param X X
//' @param y y
//' @param sigma standard error for y
//' @param betas paths of coefficients
//' @param g_index group id
//' @param grps unique groups id
//' @param K inter-group norm
//' @return a path length criterion
//' @examples 
//' \dontrun{
//' library(SA23204165)
//' library(MASS)
//' rho = 0.2
//' X <- mvrnorm(200,mu = rep(0,10), Sigma = diag(10) + rho*(1-diag(10)))
//' beta <- c(0.2,0.4,0.5,0.6,0.6,0.6,0.2,0.3,0.5,0.7)
//' g.index <- c(1,1,1,2,2,2,3,3,4,5)
//' y <- X%*%beta + rnorm(200)
//' fit <- glars(X,y,g.index)
//' criterionC(fit$X, fit$y, sd(fit$y), fit$betas, g.index, unique(g.index),c(unname(table(g.index))))
//' }
//' @export
// [[Rcpp::export]]
NumericVector criterionC(NumericMatrix X, NumericVector y,double sigma, NumericMatrix betas, NumericVector g_index, NumericVector grps, NumericVector K) {
  
    MatrixXd XX = as<MatrixXd>(X);
    VectorXd yy = as<VectorXd>(y);
    MatrixXd betas_eig = as<MatrixXd>(betas);
    VectorXd g_index_eig = as<VectorXd>(g_index);
    VectorXd grps_eig = as<VectorXd>(grps);
    VectorXd K_eig = as<VectorXd>(K);
    
    MatrixXd sigma_expand(1,betas_eig.cols()); 
    MatrixXd n_expand(1,betas_eig.cols());
    MatrixXd y_expand(yy.size(),betas_eig.cols());
    MatrixXd betas_diff(betas_eig.rows(),betas_eig.cols());
    MatrixXd group_norm(grps_eig.size(),betas_eig.cols());
    MatrixXd group_diff_norm(grps_eig.size(),betas_eig.cols());
    MatrixXd cumsum_group_diff_norm(grps_eig.size(),betas_eig.cols());
    n_expand = MatrixXd::Ones(1,betas_eig.cols()); 
    n_expand = XX.rows()*n_expand;
    for (int i = 0; i < betas.cols();i++){
      y_expand.col(i) = yy;
    }
    MatrixXd residuals =  y_expand - XX*betas_eig;
    VectorXd loss =  residuals.colwise().norm();
    MatrixXd norm_j(1,betas_eig.cols());
    for (int j = 0; j < grps_eig.size();j++){
      norm_j = MatrixXd::Zero(1,betas.cols());  
      for(int k = 0; k < betas_eig.rows();k++){
        if (grps_eig(j) == g_index_eig(k)){
          norm_j += betas_eig.row(k).cwiseAbs2();
        }
      }
      group_norm.row(j)  = norm_j.cwiseSqrt() ;
    }
    betas_diff.col(0) = betas_eig.col(0);
    for (int z = 1; z < betas_eig.cols() ;z++){
      betas_diff.col(z) = betas_eig.col(z) - betas_eig.col(z-1);
    }
    MatrixXd norm_l(1,betas_eig.cols());
    for (int l = 0; l < grps_eig.size();l++){
      norm_l = MatrixXd::Zero(1,betas_eig.cols());   
      for(int m = 0; m < betas_eig.rows();m++){
        if (grps_eig(l) == g_index_eig(m)){
          norm_l += betas_diff.row(m).cwiseAbs2();
        }
      }
      group_diff_norm.row(l)  = norm_l.cwiseSqrt() ;
    }
    cumsum_group_diff_norm.col(0) = group_diff_norm.col(0);
    for (int i = 1;i < betas_eig.cols();i++){
      cumsum_group_diff_norm.col(i) = group_diff_norm.leftCols(i+1).rowwise().sum(); 
    }
    VectorXd all = group_diff_norm.rowwise().sum();
    MatrixXd all_expand(grps_eig.size(),betas_eig.cols());
    MatrixXd K_expand(grps_eig.size(),betas_eig.cols());
    for (int j = 0; j < betas_eig.cols();j++){
      all_expand.col(j) = all;
      K_expand.col(j) = K_eig;
    }
    MatrixXd df_temp(grps_eig.size(),betas_eig.cols());
    df_temp = cumsum_group_diff_norm.cwiseQuotient(all_expand).cwiseProduct(K_expand) - cumsum_group_diff_norm.cwiseQuotient(all_expand);
    MatrixXd res1 = loss.cwiseAbs2().transpose()/pow(sigma,2) - n_expand + 2*df_temp.colwise().sum();
    NumericVector result(res1.cols());
    NumericVector count(res1.cols());
    for(int j=0;j < group_norm.cols(); j++){
      count(j) = 0;
      for(int i =0;i < group_norm.rows();i++){
        if(group_norm(i,j)!=0){
          count(j) += 2;
        }
      }
    }
    for (int k = 0; k < res1.cols();k++){
      result(k) = res1(0,k);
    }
    return result + count;
 }
