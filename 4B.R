
library(dplyr)
library(MASS)
library(ggplot2)

set.seed(1233)

local_cubic=function(DF,x,h){
  X=matrix(nrow=nrow(DF),ncol=4)
  estim_x=c()
  for (i in 1:length(x)) {
    for(p in 0:3){
        X[,(p+1)]=(DF$X-x[i])^p
    }
    W=diag(dnorm(DF$X-x[i]/h)/h)
    #estim_x[i]= (t(X)%*%W%*%X %>% ginv()) %*% (t(X) %*% W %*% DF$Y)
    estim_x[i]= t(c(1,0,0,0)) %*% (t(X)%*%W%*%X %>% ginv()) %*% (t(X) %*% W %*% DF$Y) 
  }
  
  return(estim_x)
}

n=500
X=rnorm(n,1,1)
x=seq(-4,4,l=600)

m=function(x){
  return((x-1)^2)
}

err=rnorm(n,0,0.5)

Y=m(X)+err
DF=data.frame(X,Y)


Te=local_cubic(DF,x,h)
m=m(x)

ggplot() + geom_point(data = DF,aes(x=X,y=Y),shape=1) + geom_line(aes(x,m),col="4") + 
  geom_line(aes(x,Te),col="2") + xlim(-3,4) + ylim(-2,8) + legend(c("estimato","real"),col=c("1","4"))
