## ----warning=FALSE------------------------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=4)
library(SA23204165)
library(MASS)
rho = 0.2
X <- mvrnorm(200,mu = rep(0,10), Sigma = diag(10) + rho*(1-diag(10)))
beta <- c(0.2,0.4,0.5,0.6,0.6,0.6,0.2,0.3,0.5,0.7)
g.index <- 1:length(beta)   #ordinary form
y <- X%*%beta + rnorm(200)
fit1 <- glars(X,y,g.index)
betas1 = fit1$betas
rownames(betas1) = paste('betas',1:length(beta))
colnames(betas1) = paste('path',unique(g.index))
betas1

## ----warning=FALSE------------------------------------------------------------
library(lars)  #lars package
fit2 <- lars(X,y,type = 'lar',intercept = F,normalize = T,use.Gram = F)
betas2 = fit2$beta
betas2 = betas2[-1,] #remove intercept
rownames(betas2) = paste('betas',1:length(beta))
colnames(betas2) = paste('path',unique(g.index))
betas2

## -----------------------------------------------------------------------------
cp1 = criterion(fit1$X,fit1$y,sd(fit1$y),fit1$betas, g.index, c(unname(table(g.index))))
cp2 = fit2$Cp
beta1.best = betas1[,cp1 == min(cp1)]
beta2.best = betas2[,unname(cp2[-1] == min(cp2[-1]))]
beta.ls = lm(y~scale(X)-1)$coefficient
beta.res = rbind(beta,beta.ls,beta1.best,beta2.best)
row.names(beta.res) = c('True','LS','ours','lars package')
colnames(beta.res) = colnames(beta1.best)
beta.res

## -----------------------------------------------------------------------------
path1 = fit1$path
path2 = unlist(fit2$actions)
path.res <- rbind(path1,path2)
colnames(path.res) = paste('path',unique(g.index))
rownames(path.res) = c('ours','lars package')
path.res

## -----------------------------------------------------------------------------
path(fit1)

## -----------------------------------------------------------------------------
plot(fit2)

## -----------------------------------------------------------------------------
library(microbenchmark)
comp <- microbenchmark(Cp.r =criterion(fit1$X, fit1$y, sd(fit1$y), fit1$betas, g.index, c(unname(table(g.index)))) , Cp.c = criterionC(fit1$X, fit1$y, sd(fit1$y), fit1$betas, g.index, unique(g.index),c(unname(table(g.index))))) 
comp

## -----------------------------------------------------------------------------
cp.r = criterion(fit1$X, fit1$y, sd(fit1$y), fit1$betas, g.index, c(unname(table(g.index))))
cp.c = criterionC(fit1$X, fit1$y, sd(fit1$y), fit1$betas, g.index, unique(g.index),c(unname(table(g.index))))
comp.res = rbind(cp.r,cp.c)
colnames(comp.res) = paste('path',unique(g.index))
comp.res

## -----------------------------------------------------------------------------

rho = 0.2
X <- mvrnorm(200,mu = rep(0,10), Sigma = diag(10) + rho*(1-diag(10)))
beta <- c(0,0,0,0.6,0.4,0.2,0.2,0.3,0,0)
g.index <- c(1,1,1,2,2,2,3,3,4,5)
y <- X%*%beta + rnorm(200)
fit <- glars(X,y,g.index)
path(fit)

## -----------------------------------------------------------------------------
cp = criterion(fit$X, fit$y, sd(fit$y), fit$betas, g.index, c(unname(table(g.index))))
fit$betas[,which(cp == min(cp))]

