#' @useDynLib SA23204165
#' @import MASS
#' @import stats
#' @title Calculate group version Cp criterion
#' @description Calculate model selection criterion 
#' @param X X 
#' @param y y
#' @param sigma standard error for y
#' @param betas paths of coefficients
#' @param g.index group id
#' @param K K inter-group norm 
#' @return a path length criterion
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
#' criterion(fit$X, fit$y, sd(fit$y), fit$betas, g.index, c(unname(table(g.index))))
#' }
#' @export
criterion <- function(X,y,sigma,betas,g.index,K){
  loss = apply(sweep(X%*%betas,1,y),2,my.norm)
  grps = unique(g.index)#组别
  group.norm <- function(beta){  #一列值
    group.norm.f <- function(j){ #第j组的norm
    my.norm(beta[g.index == j])
    }
  sapply(grps,group.norm.f)
  }
norms = apply(betas,2,group.norm)
df = numeric(ncol(X))
betas.diff = matrix(0,nrow(betas),ncol(betas))
betas.diff[,1] = betas[,1]
for (i in 2:ncol(betas)){
    betas.diff[,i] = betas[,i] - betas[,i-1]
  }
diff.norms = apply(betas.diff,2,group.norm)
temp = t(apply(diff.norms, 1, cumsum))
all = rowSums(diff.norms)
df = colSums(norms > 0) + colSums((temp/all)*(K-1))
loss^2/sigma^2 - nrow(X) + 2*df
}

my.norm <- function(x){
  sqrt(sum(x^2))
}