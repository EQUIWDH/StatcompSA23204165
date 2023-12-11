## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=4)
library(naniar)
library(survival)
library(MASS)
library(survminer)
library(ggpubr)
library(glmnet)

## -----------------------------------------------------------------------------
beta = c(1,0,-1,0,1,0,-1,0)
sigma <- diag(rep(1,8))    #generate variance matrix
mean <- rep(0,8)
X <- mvrnorm(200,mean,sigma)
Y <- rnorm(X %*%beta,1)
X[1:5,]

## -----------------------------------------------------------------------------
cvfit <- cv.glmnet(X,Y,family = 'gaussian',alpha = 1,nlambda = 200,intercept = FALSE)
plot(cvfit,sign.lambda = cvfit$lambda.min)

## -----------------------------------------------------------------------------
res <- cbind(coef(cvfit,cvfit$lambda.min) , append(beta,0,after = 0))
colnames(res) <- c('real','estimation')
res

## -----------------------------------------------------------------------------
data(cancer)
lung[1:5,]

## -----------------------------------------------------------------------------
# impute NAs by medians of all columns
input <- lung  %>% impute_median_if(is.numeric)
input$status[input$status == 1] = 0
input$status[input$status == 2] = 1
input[1:5,]

## -----------------------------------------------------------------------------
warnings('off')
fit <- coxph(Surv(time,status) ~.,data = input)
ggforest(fit,data = input)

## -----------------------------------------------------------------------------
women_weight <- c(88.9, 81.2, 73.3, 21.8, 63.4, 84.6, 28.4, 28.8, 28.5)
men_weight <- c(37.8, 80, 33.4, 36, 89.4, 83.3, 97.3, 81.3, 92.4)
sd(women_weight)
sd(men_weight)

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
qqnorm(women_weight)
qqline(women_weight)
qqnorm(men_weight)
qqline(men_weight)

## -----------------------------------------------------------------------------
compare_data <- data.frame( 
                group = rep(c("Woman", "Man"), each = 9),
                weight = c(women_weight,  men_weight)
                )
p = ggboxplot(compare_data,x = 'group',y ='weight',color = "group", palette = "jama",
          add = "jitter")
p + stat_compare_means()

## -----------------------------------------------------------------------------
my_sample = function(size,cdf) {
  U = runif(size)    #cdf is a matrix with variable and probability
  vars = cdf[1,]
  probs = cdf[2,]
  inds = findInterval(U,probs)
  return(vars[inds+1])
}
X <- 1:10   #constructing cdf
p <- rep(1/10,10)
cdf <- rbind(X,cumsum(p))
table(my_sample(1e5,cdf))/1e5

## -----------------------------------------------------------------------------
cdf

## -----------------------------------------------------------------------------
my.laplace = function(size) {
  U = runif(size)    
  F.inv = function(y){
    if (y >= 0 && y <= 0.5){
      return(log(2*y))
    }
    else{
      return(-log(2-2*y))
    }
  }
  y = sapply(U,F.inv)
  hist(y,probability = TRUE,main = 'Histogram of Laplace')
  y <- seq(-5,5,0.01)
  lines(y,0.5*exp(-abs(y)),lwd = 2)
}
my.laplace(500)

## -----------------------------------------------------------------------------
n <- 1e4; j<-k<-0;y <- numeric(n)
while (k < n) {
u <- runif(1)
j <- j + 1
x <- runif(1)
if (x^2 * (1-x) > 16/9*u) {
k <- k + 1
y[k] <- x
}
}
hist(y,probability = TRUE, main = 'Histogram of beta(3,2)')
y <- seq(0,1,0.01)
lines(y,12*y^2*(1-y),lwd = 2)

## -----------------------------------------------------------------------------
n <- 1e4; k<-1;y <- numeric(n)
while(k<=n){
U = runif(3,-1,1)
if (abs(U[3]) >= abs(U[2]) && abs(U[3]) >= abs(U[1])) {
  y[k] = U[2]
}
else{
  y[k] = U[3]
}
k <- k+1
}
hist(y,probability = TRUE,main = 'Histogram of fe')
y <- seq(-1,1,0.01)
lines(y,3/4*(1-y^2),lwd = 2)

## -----------------------------------------------------------------------------
set.seed(12345)
montecarlo_repeat <- function(rho,K){
d = 1 # baseline of ratio
l = rho*d
pihats = rep(0,K)
for (i in 1:K){     #repeat K times
  n <- 1e6
  X <- runif(n,0,d/2)  
  Y <- runif(n,0,pi/2)
  pihats[i] <- 2*rho/mean(l/2*sin(Y)>X)
}
return(c(mean(pihats),sd(pihats)))
}
shown = rbind(montecarlo_repeat(1,100),
              montecarlo_repeat(0.5,100)
              ,montecarlo_repeat(0.3,100))
colnames(shown) = c('mean','std')
rownames(shown) = c(1,0.5,0.3) #simulate on rho = 1,0.5,0.3
shown

## -----------------------------------------------------------------------------
1-(exp(1)^2 + 2*exp(1) - 1 - 4*(exp(1)-1)^2)/(-2*exp(1)^2 + 8*exp(1) - 6)

## -----------------------------------------------------------------------------
my_compare <- function(size){
n <- size
u <- runif(size,0,1)  
v <- runif(size,0,1)
X <- c(u,v)
Y <- exp(X)
Y1 <- exp(u)
Y2 <- exp(1-u)
var_MC <- sum((Y -  mean(Y))^2) /(2*n-1)
var_Antithetic= sum( ((Y1+Y2)/2 -  mean((Y1+Y2)/2))^2 ) /(n-1)
return(c(mean(Y),mean((Y1+Y2)/2), var_Antithetic/var_MC))
}
my_compare(1e6) # 'MC result','Antithetic result','var ratio'

## -----------------------------------------------------------------------------
library(ggplot2)
g <- function(x) return(x^2 / sqrt(2*pi)*exp(-x^2/2)*(x >= 1))
f1 <- function(x) return(exp(-x))
f2 <- function(x) return(1 / sqrt(2*pi) * exp(-x^2/2) )
t = seq(1,5,0.01)
df=data.frame(x = t,
              values = c(g(t),f1(t),f2(t)),
              func = c(rep("g(x)",length(t)), rep("f1(x)",length(t)),
                          rep("f2(x)",length(t)))
              )
  
ggplot(df, aes(x, values, col = func))+geom_line()

## -----------------------------------------------------------------------------
x1 <- rexp(1e6)  #sampling from f1
x2 <- rnorm(1e6) #sampling from f2
c(mean(g(x1)/f1(x1)), sd(g(x1)/f1(x1)), mean(g(x2)/f2(x2)), sd(g(x2)/f2(x2)))

## -----------------------------------------------------------------------------
g <- function(x) return(x^2 / sqrt(2*pi)*exp(-x^2/2)*(x >= 1))
x <- runif(1e6,1,5)
c(mean(g(x)), sd(g(x)))

## -----------------------------------------------------------------------------
g <- function(x) return(exp(-x)/(1+x^2)*(x >= 0)*(x <= 1))
stratified_est = numeric(5)
N <- 50 #replication = 50
n = 1e6
k = 5
M = n/k
estimates <- matrix(0, N, 2)
for (i in 1:N){
  # first column is stratified sample
  for (j in 1:k){
    x1 <- runif(M,(j-1)/k,j/k)
  stratified_est[j] = mean(g(x1))
  }
  estimates[i,1] = mean(stratified_est)
  # second column is importance sample 
  u <- runif(M) #f3, inverse transform method
  x2 <- -log(1 - u * (1 - exp(-1)))
  fg <- g(x2) / (exp(-x2) / (1 - exp(-1)))
  estimates[i,2] <- mean(fg)
  }
c(mean(estimates[,1]),sd(estimates[,1]),mean(estimates[,2]),sd(estimates[,2]))

## -----------------------------------------------------------------------------
alpha = c(0.01,0.05,0.1)
n = 1e3 #replication number
p <- numeric(n)
for (i in 1:n){
x <- rchisq(1e5,1) #sampling from chi^2(1)
p[i] <- t.test(x,mu = 1,alternative = "two.sided")$p.value #One sample test
}
c(mean(p <= alpha[1]),mean(p <= alpha[2]),mean(p <= alpha[3]))

## -----------------------------------------------------------------------------
alpha = c(0.01,0.05,0.1)
n = 1e3 #replication number
p <- numeric(n)
for (i in 1:n){
x <- runif(1e5,0,2) #sampling from Uniform[0,2]
p[i] <- t.test(x,mu = 1,alternative = "two.sided")$p.value #One sample test
}
c(mean(p <= alpha[1]),mean(p <= alpha[2]),mean(p <= alpha[3]))

## -----------------------------------------------------------------------------
alpha = c(0.01,0.05,0.1)
n = 1e3 #replication number
p <- numeric(n)
for (i in 1:n){
x <- rexp(1e5,1) #sampling from Exp(1)
p[i] <- t.test(x,mu = 1,alternative = "two.sided")$p.value #One sample test
}
c(mean(p <= alpha[1]),mean(p <= alpha[2]),mean(p <= alpha[3]))

## -----------------------------------------------------------------------------
alpha = 0.05
m = 1e5 #replication
mu_is_in <- numeric(n) #record if mu=2 in confidence interval
for (i in 1:m){
n = 20  #sample size
x <- rchisq(n,2)
low = mean(x) - qt(1-alpha/2,n-1) *sd(x) / sqrt(n-1) #lower bound
high = mean(x) + qt(1-alpha/2,n-1) *sd(x) /sqrt(n-1) #higher bound
mu_is_in[i] <- ifelse(low <= 2 & high >= 2,1,0) 
}
mean(mu_is_in)

## -----------------------------------------------------------------------------
library(boot)
data("aircondit")
res <- matrix(0,4,2)                       #storage 
colnames(res) = c('low','up')
row.names(res) = c('norm','basic','percent','BCa')
aircondit <-as.array(aircondit[,1])        #adapt data form
boot.MLE <- function(x,i) 1/mean(x[i])     # calculate MLE
boot.data <- boot(data = aircondit,statistic = boot.MLE,R = 1000)
ci <- boot.ci(boot.data,type=c("norm","basic","perc","bca"))
res[1,] <- c(ci$normal[2],ci$normal[3])
res[2,] <- c(ci$basic[4],ci$basic[5])
res[3,]  <- c(ci$percent[4],ci$percent[5])
res[4,] <- c(ci$bca[4],ci$bca[5])
res

## -----------------------------------------------------------------------------
library(bootstrap)
data("scor")
cor.mat <- cor(scor)  #sample correlation coef
n = dim(scor)[1]      # row length
first.princ <- max(eigen(cor.mat)$values) / sum(eigen(cor.mat)$values) #estimated first principle
jack <- numeric(dim(scor)[1])   #storage for jacknife
for (i in 1:n){
  mat.jack <- cor(scor[-i,])      
  jack[i] <- max(eigen(mat.jack)$values) / sum(eigen(mat.jack)$values) #jacknife estimation
}
bias.jack <- (n-1)*(mean(jack)-first.princ)
se.jack <- sqrt((n-1)*mean((jack-first.princ)^2))
round(c(original= first.princ, bias.jack=bias.jack, se.jack = se.jack),4)


## -----------------------------------------------------------------------------
library(DAAG)
data(ironslag)
magnetic <- ironslag$magnetic
chemical <-  ironslag$chemical
n = dim(ironslag)[1]
e1.one <- e2.one <- e3.one <- e4.one <- numeric(n)
e1.two <- e2.two <- e3.two <- e4.two <-  numeric(n*(n-1)/2)

# leave-one-out
for (k in 1:n) {
y <- magnetic[-k]
x <- chemical[-k]
J1 <- lm(y ~ x)
yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
e1.one[k] <- magnetic[k] - yhat1
J2 <- lm(y ~ x + I(x^2))
yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] +
J2$coef[3] * chemical[k]^2
e2.one[k] <- magnetic[k] - yhat2
J3 <- lm(log(y) ~ x)
logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
yhat3 <- exp(logyhat3)
e3.one[k] <- magnetic[k] - yhat3
J4 <- lm(log(y) ~ log(x))
logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[k])
yhat4 <- exp(logyhat4)
e4.one[k] <- magnetic[k] - yhat4
}

#leave-two-out
k = 1
for (i in 1:(n-1)){
  for (j in (i+1):n){
    y <- magnetic[-c(i,j)]
    x <- chemical[-c(i,j)]
    
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(i,j)]
    e1.two[k] <- sum((magnetic[c(i,j)] - yhat1)^2)
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(i,j)] +
    J2$coef[3] * chemical[c(i,j)]^2
    e2.two[k] <- sum((magnetic[c(i,j)] - yhat2)^2)
    
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(i,j)]
    yhat3 <- exp(logyhat3)
    e3.two[k] <- sum((magnetic[c(i,j)] - yhat3)^2)
    
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(i,j)])
    yhat4 <- exp(logyhat4)
    e4.two[k] <- sum((magnetic[c(i,j)] - yhat4)^2)
    
    k <- k + 1
  }
}
res <- rbind(c(mean(e1.one^2), mean(e2.one^2), mean(e3.one^2), mean(e4.one^2)),
                c(sum(e1.two)/(n*(n-1)), sum(e2.two)/(n*(n-1)), sum(e3.two)/(n*(n-1)), sum(e4.two)/(n*(n-1))))
colnames(res) <- c('e1','e2','e3','e4')
row.names(res) <- c('leave-one-out','leave-two-out')
res

## -----------------------------------------------------------------------------
cvm <- function(x,y){      #calculate cramer-von Mises statistics
  n = length(x)
  m = length(y)
  f.n <- ecdf(x)
  g.m <- ecdf(y)
  return(n*m / (n+m)^2 * (sum((f.n(x) - g.m(x))^2) + sum((f.n(y) - g.m(y))^2) ))
}
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

R <- 10000 #number of replication
z <- c(x, y)    #pooling 
K <- length(z)
perm <- numeric(R) #storage for replicates
options(warn = -1)
test0 <- cvm(x,y)
for (i in 1:R) {
#generate indices k for the first sample
k <- sample(K, size = 14, replace = FALSE)
x.perm <- z[k]
y.perm <- z[-k] 
perm[i] <- cvm(x.perm,y.perm)
}
p <- mean(c(test0, perm ) >= test0)
options(warn = 0)
p

## -----------------------------------------------------------------------------
count5test <- function(x, y) {  #function to do test
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}

n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0  # mean
sigma1 <- sigma2 <- 1 # variance
x <- rnorm(n1,mu1,sigma1)
y <- rnorm(n2,mu2,sigma2)
z <- c(x,y)
m <- 9999  #replication
perm <- numeric(m) #storage for perm
for (i in 1:m) {
#generate indices k for the first sample
k <- sample(n1+n2, size = n1, replace = FALSE)
x.perm <- z[k]
y.perm <- z[-k] 
perm[i] <- count5test(x.perm,y.perm)
}
mean(perm)

## -----------------------------------------------------------------------------
answer <- function(N,b1,b2,b3,f0){
  logistics <- function(a){
  return(mean(exp(b1*x1+b2*x2+b3*x3+a) / (1 + exp(b1*x1+b2*x2+b3*x3+a))) - f0)
}
    x1 <- rpois(N,1)
    x2 <- rexp(N)
    x3 <- rbinom(N,1,0.5)
    return(unlist(uniroot(logistics,c(-30,30))))
}
f <- c(0.1, 0.01, 0.001, 0.0001)
ans <- matrix(0,4,5)
for (i in 1:4){
ans[i,] <- answer(10e6,0,1,-1,f[i])
}
rownames(ans) <- f
colnames(ans) <- names(answer(10e6,0,1,-1,f[i]))
ans

## -----------------------------------------------------------------------------
plot(-log(f),ans[,1],xlab = '-log(f0)',ylab = 'alpha',main = 'result')

## -----------------------------------------------------------------------------
sigmas <- c(0.1,1,4,16)
laplace <- function(x){
  return(1/2*exp(-abs(x)))
}
my.chain <- function(sigma){
  N <- 5000
  xt  <- numeric(N)
  accept <- 0
  xt[1] <- rnorm(1,0,0.1) 
  for (i in 2:N){
    y <- rnorm(1,xt[i-1],sigma) # sample from Proposal distribution
    u <- runif(1)
    if (u <= laplace(y) / laplace(xt[i-1])){
      accept <- accept + 1
      xt[i] <- y
      
    }
    else{
      xt[i] <- xt[i-1]
    }
  }
  return(list(x = xt, ac = accept/N))
}

accepts <- numeric(4)
names(accepts) <- sprintf("sigma = %.2f", sigmas) #storage for accepts and quantile for sigmas
quant.table <- matrix(0,9,5)
colnames(quant.table) <- c('theory',sprintf("sigma = %.2f", sigmas))
row.names(quant.table) <- paste0(c(1:9)*10, "%")
quant.table[,1] <-VGAM::qlaplace(c(1:9)/10)
i <- 1
for (sigma in sigmas){
  tmp <- my.chain(sigma)
  accepts[i] <- tmp$ac
  quant.table[,i+1] <- quantile(tmp$x,c(1:9)/10)
  i <- i+1
}
quant.table

## -----------------------------------------------------------------------------
accepts

## -----------------------------------------------------------------------------
N <- 5000
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
rho = 0.9
gibbs <- matrix(0,5000,2)
colnames(gibbs) <- c('xt','yt')
gibbs[1,] <- c(0,0)
for (i in 2:5000){
  gibbs[i,1] <- rnorm(1,mu1 + rho*sigma1/sigma2*(gibbs[i-1,2] - mu2), sqrt(1 - rho^2)*sigma1)
  gibbs[i,2] <- rnorm(1,mu2 + rho*sigma2/sigma1*(gibbs[i,1] - mu1), sqrt(1 - rho^2)*sigma2)
}
par(mfrow = c(1,2))
plot(gibbs[,1][501:5000],type = 'l',xlab = 'iter', main = 'xt',
  ylab = "value", ylim = range(gibbs[,1][501:5000]))
plot(gibbs[,2][501:5000],type = 'l',xlab = 'iter',main = 'yt',
  ylab = "value", ylim = range(gibbs[,1][501:5000]))


## -----------------------------------------------------------------------------
reg <- lm(gibbs[,1]~gibbs[,2])
hist(residuals(reg),probability = TRUE, xlab = 'residuals',main = 'residuals of regression')

## -----------------------------------------------------------------------------
library(coda)

Rayleigh <- function(x, sigma) {
if (any(x < 0)) return (0)
stopifnot(sigma > 0)
return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

Gelman.Rubin <- function(ps) {
# ps[i,j] is the statistic ps(X[i,1:j])
# for chain in i-th row of X
ps <- as.matrix(ps)
n <- ncol(ps)
k <- nrow(ps)
ps.means <- rowMeans(ps) #row means
B <- n * var(ps.means) #between variance est.
ps.w <- apply(ps, 1, "var") #within variances
W <- mean(ps.w) #within est.
v.h <- W*(n-1)/n + (B/n) #upper variance est.
r.h <- v.h / W #G-R statistic
return(r.h)
}

my.chains <- function(K,N,sigma){ #input is length of chain and sigma 
  
sigma <- sigma
x <- matrix(0,K,N)
gb.stat <- numeric(N)
x[,1] <- rchisq(K, df=1)
for (j in 2:N) {
  u <- runif(1)
  for (i in 1:K){
    y <- rchisq(1, df = x[i,j-1])
    num <- Rayleigh(y, sigma) * dchisq(x[i,j-1], df = y)
    den <- Rayleigh(x[i,j-1], sigma) * dchisq(y, df = x[i,j-1])
    if (u <= num/den) {
     x[i,j] <- y 
    }
   else{
   x[i,j] <- x[i,j-1]
   }
}
  
  ps <- t(apply(x[,1:j], 1 , cumsum)) 
  for (i in 1:nrow(ps)){
    ps[i,] <- ps[i,] / (1:j)
  }
  gb.stat[j-1]<-Gelman.Rubin(ps)
    
}
 return(list(x = x,gb = gb.stat))
 
}


res <- my.chains(5,10000,0.5)
chain1 <- as.mcmc.list(as.mcmc(res$x[1,]))
chain2 <- as.mcmc.list(as.mcmc(res$x[2,]))
chain3 <- as.mcmc.list(as.mcmc(res$x[3,]))
chain4 <- as.mcmc.list(as.mcmc(res$x[4,]))
chain5 <- as.mcmc.list(as.mcmc(res$x[5,]))

all.chains <- rbind(chain1,chain2,chain3,chain4,chain5)

gelman.plot(all.chains)


## -----------------------------------------------------------------------------
gelman.diag(all.chains)

## -----------------------------------------------------------------------------
u = c(11,8,27,13,16,0,23,10,24,2) #left points
v = c(12,9,28,14,17,1,24,11,25,3) #right points
eps <- 1e8
n = length(u)

g.MLE <- function(lambda){
  return(sum((u*exp(-lambda*u) - v*exp(-lambda*v))/(exp(-lambda*u) - exp(-lambda*v)))  )
}
iter.func <- function(lambda){
  return(1 / mean( (u*exp(-lambda*u) - v*exp(-lambda*v))/(exp(-lambda*u) - exp(-lambda*v))  ))
}
lambda0 <- 5
lambda <- iter.func(lambda0)
diff <- abs(lambda - lambda0)
while(diff > eps){
  lambda <- iter.func(lambda0)
  diff <- abs(lambda0 - lambda)
  lambda0 <- lambda
}
c(uniroot(g.MLE,c(0,1))$root,lambda)  #MLE and EM result


## -----------------------------------------------------------------------------
library(boot)
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)

my.morra <- function(mat) {
#solve the two player zero-sum game by simplex method
#optimize for player 1, then player 2
#maximize v subject to ...
#let x strategies 1:m, and put v as extra variable
#A1, the <= constraints
#
A <- (mat- min(mat)) / max(mat)
m <- nrow(A)
n <- ncol(A)
iter.times <- n^3
a <- c(rep(0, m), 1) #objective function
A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
b3 <- 1
sx <- boot::simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,maxi=TRUE, n.iter=iter.times)
a <- c(rep(0, n), 1) #objective function
A1 <- cbind(A, rep(-1, m)) #constraints <=
b1 <- rep(0, m)
A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
b3 <- 1
sy <- boot::simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=FALSE, n.iter=iter.times)
solution <- list("mat" = A * max(mat) + min(mat),
"x" = sx$soln[1:m],
"y" = sy$soln[1:n],
"v" = sx$soln[m+1] * max(mat) + min(mat))
solution
}

B <- A + 2
s <- my.morra(B)
round(cbind(s$x, s$y), 7)

## -----------------------------------------------------------------------------
test <- list(x = 1,y = c(1,2))
a = as.vector(test)
a # it seems not working, we must use unlist

## -----------------------------------------------------------------------------
a <- c(1, 2, 3, 4)
dim(a)  #return NULL because vector is a 0-array

## -----------------------------------------------------------------------------
b <- matrix(1:9,3,3)
is.matrix(b)
is.array(b)  #matrix is a speciall array with only 2 dims

## -----------------------------------------------------------------------------
my.df <- data.frame(col1 =rep(TRUE,5),col2 = 1:5, col3 = rnorm(5), col4 = rep('str',5))
as.matrix(my.df)  

## -----------------------------------------------------------------------------
a = data.frame(x = c(1,2,3,4))
a.T = t(a)
a[,-1]  #Actuall, we can get such dataframes

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------



## -----------------------------------------------------------------------------
num.df<-data.frame(matrix(sample(2:10,30,replace = TRUE),5,6))
vapply(num.df,sd,FUN.VALUE = c(se = 0))

## -----------------------------------------------------------------------------
mix.df<-data.frame(x = rnorm(5), y = sample(c('a','b','c'),5,replace = TRUE),z = runif(5))
is.num <- function(x){    #input is mix.df
  vapply(x, is.numeric, FUN.VALUE = TRUE)
}
vapply(mix.df[,is.num(mix.df)], mean , FUN.VALUE = 0)

## -----------------------------------------------------------------------------
library(Rcpp)
gibbs.r <- function(a,b,n,N){
  output<-matrix(0,N,2)
  for (i in 2:N){
    output[i,1] <- rbinom(1,n,output[i-1,2])
    output[i,2] <- rbeta(1,output[i,1]+a,n-output[i,1]+b)
  }
  output
}

sourceCpp('./gibbs.cpp')
res <- cbind(gibbs.r(3,4,10,1000), gibbs(3,4,10,1000))
colnames(res) <- c('gibbs.r x', 'gibbs.r y', 'gibbs.c x', 'gibbs.c y')
round(res,3)[1:20,]

## -----------------------------------------------------------------------------
library(microbenchmark)
comp <- microbenchmark(gibbs.r=gibbs.r(3,4,10,1000), gibbs.c =  gibbs(3,4,10,1000))
summary(comp)

