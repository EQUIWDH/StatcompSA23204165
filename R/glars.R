#' @useDynLib SA23204165
#' @import MASS
#' @importFrom pracma gramSchmidt
#' @title  A function for constructing path of group LARS
#' @description A method for constructing path of group LARS
#' @param X The observed patients number  
#' @param y The  features dimension of covariates
#' @param g.index The pairwise correlation in design matrix
#' @param K Inter-group norm, usually take p_jI_j.
#' @param normlized The pairwise correlation type of features
#' @param gram The true coefficient provided
#' @return a list with generated coefficient path and add path 
#' @examples
#' \dontrun{
#' library(SA23204165)
#' library(MASS)
#' rho = 0.2
#' X <- mvrnorm(200,mu = rep(0,10), Sigma = diag(10) + rho*(1-diag(10)))
#' beta <- c(0.2,0.4,0.5,0.6,0.6,0.6,0.2,0.3,0.5,0.7)
#' g.index <- c(1,1,1,2,2,2,3,3,4,5)
#' y <- X%*%beta + rnorm(200)
#' fit <- glars(X,y,g.index)
#' fit$betas
#' }
#' @export
glars <- function(X,y,g.index,K = c(unname(table(g.index))),normlized = TRUE,gram = TRUE){  
  if (length(g.index) == ncol(X)){                               
    colnames(X) = g.index                                          
    gamma <- numeric(length(g.index))
    grps = unique(g.index)   
    A = logical(length(grps)) 
    add.order <- numeric(length(grps)) 
    k = 1 
    alpha <- beta <- 0
    beta.path <- matrix(0,ncol(X),length(grps)) 
    r <-y   
    if (gram){
      for (col in grps){
        if (sum(g.index == col) > 1){
          X[,g.index == col] = gramSchmidt(X[,g.index == col])$Q
        }
      }
    }
    if (normlized){
      non.zeros.cols = apply(X!=0,2,sum) != 0
      X[,non.zeros.cols] = scale(X[,non.zeros.cols])
      X[is.na(X)] = 0 
    }
    corr0 = as.vector(t(X) %*% r)
    corr.k<- function(k){  
      sum(corr0[g.index == k]^2/K[which(grps == k)])
    }
    corr = sapply(grps,corr.k)
    A[corr == max(corr)] = TRUE 
    add.order[corr == max(corr)] = 1 
    while (alpha != 1){
      col = logical(length(g.index))
      for (gr in grps[A]){
        col = col | g.index == gr
      }
      gamma[col] = ginv(t(X[,col])%*%X[,col])%*%t(X[,col])%*%r
      if (length(grps[A]) > 1){
        j = sample(grps[A],1)  
      }else{
        j = grps[A]
      }
      ajc <- numeric(length(grps[!A]))  
      if (sum(A) != length(grps)){
        for (i in 1:length(grps[!A])){
          jc = grps[!A][i]
          solve.a <- function(alpha.j){   
            sum((t(X[,g.index == j]) %*%(r - alpha.j*X%*%gamma))^2)/K[which(grps == j)] - sum((t(X[,g.index == jc]) %*%(r - alpha.j*X%*%gamma))^2)/K[which(grps == jc)] 
          }
          if(solve.a(0)*solve.a(1) <0){
            ajc[i] <- uniroot(solve.a,c(0,1),tol = 1e-6)$root
          }
          if(isTRUE(all.equal(solve.a(0),0))){
            ajc[i] = 0
          }
          if(isTRUE(all.equal(solve.a(1),0))){
            ajc[i] = 1
          }
        }
        alpha = min(ajc)
        add = grps[!A][ajc == min(ajc)]  
        if (length(add) > 2 & alpha == 1){ 
          beta = beta + alpha*gamma
          beta.path[,k] = beta
          r = y - X%*%beta
          k = k+1
          beta.path = beta.path[,colSums(beta.path^2) !=0]
          grps[order(add.order)]
          return(list(X = X,y = y,index = g.index,K = K,betas = beta.path,path = grps[order(add.order)]))
        }
        A[grps == add] = TRUE
        add.order[grps == add] = k+1
        beta = beta + alpha*gamma
        beta.path[,k] = beta
        r = y - X%*%beta
        k = k+1
        
      }
      else{
        alpha = 1
        beta.path[,k] = beta + alpha*gamma
      }
      
    }
    
  }
  
  list(X = X,y = y,index = g.index,K = K,betas = beta.path[,colSums(beta.path^2) !=0],path = grps[order(add.order)][colSums(beta.path^2) !=0])
  
}