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
