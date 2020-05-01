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

## mcmc
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

## importance sampling
sum_numerator<-rep(0,5)
sum_denominator<-0
cc<-coef(glm.out)
ee<-sweep(matrix(rnorm(5*100000),100000,5),2,cc,"+")
for(iter in 1:100000){
  ww<-exp(logl(ee[iter,])+sum((ee[iter,]-cc)^2)/2)
  sum_numerator<-sum_numerator+ee[iter,]*ww
  sum_denominator<-sum_denominator+ww
}
sum_numerator/sum_denominator

## mcmc package
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

### Logistic mixed model

```
## mixed model logistic regression

set.seed(10)
nind<-100
nobs<-200
beta0_true<-1
beta1_true<-2
alpha_true<-rnorm(nind)
xx<-matrix(rnorm(nind*nobs),nind,nobs)
yy<-matrix(rbinom(nind*nobs,size=1,prob=plogis(beta0_true+beta1_true*xx+alpha_true)),nind,nobs)

###########################
# variance of alpha known #
###########################

niter<-5000
all_beta0<-rep(NA,niter)
all_beta1<-rep(NA,niter)
all_alpha<-matrix(NA,nind,niter)
beta0_est<-1
beta1_est<-2
alpha_est<-alpha_true

logpy_full<-function(beta0,beta1){
  eta<-c(beta0+beta1*xx+alpha_est)
  logp<-ifelse(eta<0,eta-log1p(exp(eta)),-log1p(exp(-eta)))
  logq<-ifelse(eta<0,-log1p(exp(eta)),-eta-log1p(exp(-eta)))
  logl<-sum(logp[c(yy)==1])+sum(logq[c(yy)==0])
  return(logl-beta0^2/32-beta1^2/32)
}

logpyi_full<-function(idx,alphai){
  etai<-beta0_est+beta1_est*xx[idx,]+alphai
  logpi<-ifelse(etai<0,etai-log1p(exp(etai)),-log1p(exp(-etai)))
  logqi<-ifelse(etai<0,-log1p(exp(etai)),-etai-log1p(exp(-etai)))
  yyi<-yy[idx,]
  logli<-sum(logpi[yyi==1])+sum(logqi[yyi==0])
  return(logli-alphai^2/2)
}

for(iter in 1:niter){
  cat(iter,"-")
  beta0_prop<-rnorm(1,beta0_est,0.2)
  beta1_prop<-rnorm(1,beta1_est,0.2)
  HR<-exp(logpy_full(beta0_prop,beta1_prop)-logpy_full(beta0_est,beta1_est))
  if(runif(1)<HR){
    beta0_est<-beta0_prop
    beta1_est<-beta1_prop
  }
  all_beta0[iter]<-beta0_est
  all_beta1[iter]<-beta1_est
  for(ii in 1:nind){
    alpha_prop<-rnorm(1,alpha_est[ii],0.2)
    HR<-exp(logpyi_full(ii,alpha_prop)-logpyi_full(ii,alpha_est[ii]))
    if(runif(1)<HR){
      alpha_est[ii]<-alpha_prop
    }
    all_alpha[ii,iter]<-alpha_est[ii]
  }
}

plot(all_beta0,type="l")
plot(all_alpha[1,1:100],type="l")

mean(all_beta0[3001:niter])
mean(all_beta1[3001:niter])

plot(all_beta0)

var(rowMeans(all_alpha[,3001:niter]))
mean(rowMeans(all_alpha[,3001:niter]))

rowMeans(all_alpha[,501:niter])-alpha_true
plot(rowMeans(all_alpha[,501:niter]))
points(alpha_true)
```

```
## mixed model logistic regression

set.seed(10)
nind<-100
nobs<-200
beta0_true<-1
beta1_true<-2
sigma2_true<-4
alpha_true<-rnorm(nind,0,sqrt(sigma2_true))
xx<-matrix(rnorm(nind*nobs),nind,nobs)
yy<-matrix(rbinom(nind*nobs,size=1,prob=plogis(beta0_true+beta1_true*xx+alpha_true)),nind,nobs)

#############################
# variance of alpha unknown #
#############################

niter<-5000
all_beta0<-rep(NA,niter)
all_beta1<-rep(NA,niter)
all_sigma2<-rep(NA,niter)
all_alpha<-matrix(NA,nind,niter)
beta0_est<-1
beta1_est<-2
sigma2_est<-1
alpha_est<-alpha_true

logpy_full<-function(beta0,beta1){
  eta<-c(beta0+beta1*xx+alpha_est)
  logp<-ifelse(eta<0,eta-log1p(exp(eta)),-log1p(exp(-eta)))
  logq<-ifelse(eta<0,-log1p(exp(eta)),-eta-log1p(exp(-eta)))
  logl<-sum(logp[c(yy)==1])+sum(logq[c(yy)==0])
  return(logl-beta0^2/32-beta1^2/32)
}

logpyi_full<-function(idx,alphai){
  etai<-beta0_est+beta1_est*xx[idx,]+alphai
  logpi<-ifelse(etai<0,etai-log1p(exp(etai)),-log1p(exp(-etai)))
  logqi<-ifelse(etai<0,-log1p(exp(etai)),-etai-log1p(exp(-etai)))
  yyi<-yy[idx,]
  logli<-sum(logpi[yyi==1])+sum(logqi[yyi==0])
  return(logli-alphai^2/2/sigma2_est)
}

for(iter in 1:niter){
  cat(iter,"-")
  beta0_prop<-rnorm(1,beta0_est,0.1)
  beta1_prop<-rnorm(1,beta1_est,0.1)
  HR<-exp(logpy_full(beta0_prop,beta1_prop)-logpy_full(beta0_est,beta1_est))
  if(runif(1)<HR){
    beta0_est<-beta0_prop
    beta1_est<-beta1_prop
  }
  all_beta0[iter]<-beta0_est
  all_beta1[iter]<-beta1_est
  for(ii in 1:nind){
    alpha_prop<-rnorm(1,alpha_est[ii],0.2)
    HR<-exp(logpyi_full(ii,alpha_prop)-logpyi_full(ii,alpha_est[ii]))
    if(runif(1)<HR){
      alpha_est[ii]<-alpha_prop
    }
    all_alpha[ii,iter]<-alpha_est[ii]
  }
  sigma2_est<-1/rgamma(1,1+nind/2,1+sum(alpha_est^2)/2)
  all_sigma2[iter]<-sigma2_est
}

plot(all_beta0[2001:niter],type="l")
plot(all_alpha[1,2001:niter],type="l")
plot(all_sigma2[2001:niter],type="l")

mean(all_sigma2)

mean(all_beta0[3001:niter])
mean(all_beta1[3001:niter])
```

Refenreces:
[1] Statistical Ecology - Ruth King - Handbook of Markov Chain Monte Carlo
[2] Sampling from the posterior distribution in generalized linear mixed models - Dani Gamerman


