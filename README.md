# learn_Bayesian_MCMC

```

# Metropolis-Hastings Update
# current x, next is q(x,.)
# Hastings ratio: r(x,y)=[h(y)*q(y,x)]/[h(x)*q(x,y)]
# Accept with prob: a(x,y)=min(1,r(x,y))

niter<-1000
xx_vec<-rep(NA,niter)
hh<-function(xx)exp(-(xx-10)^2/2)
xx<-0
for(iter in 1:niter){
  yy<-rnorm(1,xx)
  # HR<-(hh(yy)*exp(-(yy-xx)^2/2))/(hh(xx)*exp(-(yy-xx)^2/2))
  HR<-hh(yy)/hh(xx)
  if(runif(1)<HR){
    xx<-yy
  }
  xx_vec[iter]<-xx
}

hist(xx_vec)
mean(xx_vec)
```

```
library(mcmc)
data(logit)
glm.out <- glm(y ~ x1 + x2 + x3 + x4, data = logit,
           family = binomial(), x = TRUE)
summary(glm.out)

x <- glm.out$x
y <- glm.out$y

logl <- function(beta) {
  eta <- as.numeric(x %*% beta)
  logp <- ifelse(eta < 0, eta - log1p(exp(eta)),- log1p(exp(- eta)))
  logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
  logl <- sum(logp[y == 1]) + sum(logq[y == 0])
  return(logl - sum(beta^2)/8)
}

niter<-100000
all_beta_est<-matrix(NA,niter,5)
beta_est<-rep(0,5)
for(iter in 1:niter){
  beta_prop<-beta_est+rnorm(5)/10
  HR<-exp(logl(beta_prop)-logl(beta_est))
  if(runif(1)<HR){
    beta_est<-beta_prop
  }
  all_beta_est[iter,]<-beta_est
}
colMeans(all_beta_est[50001:niter,])

coef(glm.out)

lupost <- function(beta, x, y, ...) {
  eta <- as.numeric(x %*% beta)
  logp <- ifelse(eta < 0, eta - log1p(exp(eta)),- log1p(exp(- eta)))
  logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
  logl <- sum(logp[y == 1]) + sum(logq[y == 0])
  return(logl - sum(beta^2) / 8)
}
beta.init <- as.numeric(coefficients(glm.out))
out <- metrop(lupost, beta.init, nbatch=5000, blen=1000, scale=0.2, x = x, y = y)
out$accept
dim(out$batch)
colMeans(out$batch[1001:5000,])
```
