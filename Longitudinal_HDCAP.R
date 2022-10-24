#################################################
# High-dimensional longitudinal Covariance regression
# linear shrinkage on the covariance matrix
#################################################

library("MASS")       # general inverse of a matrix
# library("nloptr")     # non-linear optimization
library("multigroup") # common PCA

# library("glmnet")     # lasso package

# library("lme4")
library("nlme")
#################################################

###########################################
# Bias corrected confidence interval for bootstrap
BC.CI <- function(theta,sims,conf.level=0.95) 
{
  low <- (1 - conf.level)/2
  high <- 1 - low
  z.inv <- length(theta[theta < mean(theta)])/sims
  z <- qnorm(z.inv)
  U <- (sims - 1) * (mean(theta) - theta)
  top <- sum(U^3)
  under <- (1/6) * (sum(U^2))^{3/2}
  a <- top/under
  lower.inv <- pnorm(z + (z + qnorm(low))/(1 - a * (z + qnorm(low))))
  lower2 <- lower <- quantile(theta, lower.inv)
  upper.inv <- pnorm(z + (z + qnorm(high))/(1 - a * (z + qnorm(high))))
  upper2 <- upper <- quantile(theta, upper.inv)
  return(c(lower, upper))
}
###########################################

#################################################
# standardized Frobenius norm
norm.F.std<-function(A1,A2=NULL)
{
  p<-nrow(A1)
  
  if(is.null(A2))
  {
    return(sqrt(sum(diag(A1%*%t(A1)))/p))
  }else
  {
    return(sum(diag(A1%*%t(A2)))/p)
  }
}

# LW linear shrinkage of covariance matrix
cov.ls<-function(Y)
{
  # Y: data
  
  n<-nrow(Y)
  p<-ncol(Y)
  
  # demean of Y
  Y<-scale(Y,center=TRUE,scale=FALSE)
  
  # sample covariance matrix
  S<-cov(Y)*(n-1)/n
  
  Ip<-diag(rep(1,p))
  
  m<-norm.F.std(S,diag(rep(1,p)))
  d2<-(norm.F.std(S-m*diag(rep(1,p))))^2
  
  b2.bar<-mean(apply(Y,1,function(x){return((norm.F.std(x%*%t(x)-S))^2)}))/n
  
  b2<-min(b2.bar,d2)
  
  a2<-d2-b2
  
  return(b2*m*Ip/d2+a2*S/d2)
}
#################################################

#################################################
# objective function
obj.func<-function(X,nT,Sigma,gamma,beta0.vec,beta0,sigma2,beta1)
{
  # X: covariates (longitudinal)
  # nT: # of observations of each subject at each visit
  # Sigma: an estimate of the covariance matrix
  # gamma: linear projection
  # beta0.vec: n by 1 vector, the random intercept matrix
  # beta0: mean of beta0.vec
  # sigma2: variance of beta0.vec
  # beta1: (q-1) by 1 vector of beta1
  
  n<-length(X)
  nV<-sapply(X,nrow)
  
  ll1<-0
  ll2<-0
  for(i in 1:n)
  {
    for(v in 1:nV[i])
    {
      u1<-(t(X[[i]][v,])%*%c(beta0.vec[i],beta1))[1,1]
      ll1<-ll1+(u1+(t(gamma)%*%Sigma[[i]][,,v]%*%gamma)[1,1]*exp(-u1))*nT[i,v]/2
    }
    ll2<-ll2+(log(sigma2)/2+(beta0.vec[i]-beta0)^2/(2*sigma2))
  }
  
  return(ll1+ll2)
}

# Linear shrinkage estimator of covariance matrix
# shrinkage parameter constant over subjects and visits
cov.ls.const<-function(Y,X,gamma,beta0.vec,beta1)
{
  # Y: outcome list (longitudinal)
  # X: covariates (longitudinal)
  # gamma: linear projection
  # beta0.vec: n by 1 vector, the random intercept matrix
  # beta1: (q-1) by 1 vector of beta1
  
  n<-length(Y)
  nV<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])
  
  # estimate of mu parameter
  # sample covariance matrix
  mu<-0
  nT<-matrix(NA,n,max(nV))
  S.iv<-vector("list",length=n)
  for(i in 1:n)
  {
    S.iv[[i]]<-array(NA,c(p,p,nV[i]))
    for(v in 1:nV[i])
    {
      nT[i,v]<-nrow(Y[[i]][[v]])
      # sample covariance matrix
      S.iv[[i]][,,v]<-cov(Y[[i]][[v]])*(nT[i,v]-1)/nT[i,v]
      
      u1<-(t(X[[i]][v,])%*%c(beta0.vec[i],beta1))[1,1]
      mu<-mu+exp(u1)
    }
    mu<-mu/nV[i]
  }
  mu<-mu/(n*(t(gamma)%*%gamma)[1])
  
  # estimate of delta, psi and phi
  hat.delta.iv=hat.psi.iv=hat.phi.iv<-matrix(NA,n,max(nV))
  hat.delta=hat.psi=hat.phi<-0
  for(i in 1:n)
  {
    for(v in 1:nV[i])
    {
      score.tmp<-(t(gamma)%*%S.iv[[i]][,,v]%*%gamma)[1,1]
      u1<-(t(X[[i]][v,])%*%c(beta0.vec[i],beta1))[1,1]
      # delta
      hat.delta.iv[i,v]<-(score.tmp-mu*(t(gamma)%*%gamma)[1,1])^2
      # psi
      hat.psi.iv[i,v]<-min((score.tmp-exp(u1))^2/nT[i,v],hat.delta.iv[i,v])
      # hat.psi.iv[i,v]<-min((score.tmp-exp(u1))^2,hat.delta.iv[i,v])
      # phi
      hat.phi.iv[i,v]<-hat.delta.iv[i,v]-hat.psi.iv[i,v]
      
      hat.delta<-hat.delta+hat.delta.iv[i,v]
      hat.psi<-hat.psi+hat.psi.iv[i,v]
      hat.phi<-hat.phi+hat.phi.iv[i,v]
    }
    hat.delta<-hat.delta/nV[i]
    hat.psi<-hat.psi/nV[i]
    hat.phi<-hat.phi/nV[i]
  }
  hat.delta<-hat.delta/n
  hat.psi<-hat.psi/n
  hat.phi<-hat.phi/n
  
  # shrinkage estimator of covariance matrix
  S.ls<-vector("list",length=n)
  for(i in 1:n)
  {
    S.ls[[i]]<-array(NA,c(p,p,nV[i]))
    for(v in 1:nV[i])
    {
      S.ls[[i]][,,v]<-(hat.psi*mu/hat.delta)*diag(rep(1,p))+(hat.phi/hat.delta)*S.iv[[i]][,,v]
    }
  }
  
  re<-list(S=S.ls,mu=mu,delta=hat.delta,psi=hat.psi,phi=hat.phi,delta.mat=hat.delta.iv,psi.mat=hat.psi.iv,phi.mat=hat.phi.iv,rho1=hat.psi*mu/hat.delta,rho2=hat.phi/hat.delta)
  
  return(re)
}
#################################################

#################################################
# given gamma, estimate beta
cap.cov_beta<-function(Y.cov,X,Y.n,gamma,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE)
{
  # Y.cov: a list of covariance matrix of Y
  # X: covariates (longitudinal)
  # Y.n: n by nV matrix of sample size
  # gamma: linear projection
  
  #--------------------------------------------
  n<-length(Y.cov)
  nV<-sapply(X,nrow)
  p<-dim(Y.cov[[1]])[1]
  
  q<-ncol(X[[1]])
  
  nT<-Y.n
  for(i in 1:n)
  {
    if(is.null(colnames(X[[i]]))==TRUE)
    {
      colnames(X[[i]])<-paste0("X",0:(q-1))
    }
  }
  #--------------------------------------------
  
  #--------------------------------------------
  # Initial value: estimate covariance matrix for each subject
  S.iv<-Y.cov
  score<-matrix(NA,n,max(nV))
  for(i in 1:n)
  {
    for(v in 1:nV[i])
    {
      score[i,v]<-t(gamma)%*%S.iv[[i]][,,v]%*%gamma
    }
  }
  #--------------------------------------------
  
  #--------------------------------------------
  # inintial estimate using a mixed effects model
  dtmp<-NULL
  for(i in 1:n)
  {
    dtmp<-rbind(dtmp,data.frame(ID=rep(i,nV[i]),score=log(score[i,1:nV[i]]),X[[i]][,-1]))
  }
  colnames(dtmp)<-c("ID","score",colnames(X[[1]])[-1])
  eval(parse(text=paste0("fit.tmp<-lme(score~",paste(colnames(X[[1]])[-1],collapse="+"),",random=~1|ID,data=dtmp,control=lmeControl(opt='optim'))")))
  coef.fix<-fit.tmp$coefficients$fixed
  beta1.new<-as.matrix(coef.fix[-1])
  beta0.vec.new<-coef.fix[1]+fit.tmp$coefficients$random$ID
  beta0.new<-mean(beta0.vec.new)
  sigma2.new<-mean((beta0.vec.new-beta0.new)^2)
  #--------------------------------------------
  
  if(trace)
  {
    beta1.trace<-cbind(beta1.new)
    beta0.vec.trace<-cbind(beta0.vec.new)
    beta0.trace<-cbind(beta0.new)
    sigma2.trace<-cbind(sigma2.new)
    
    obj<-obj.func(X,nT,S.iv,gamma,beta0.vec.new,beta0.new,sigma2.new,beta1.new)
  }
  
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update random beta0
    beta0.vec.upd<-rep(NA,n)
    for(i in 1:n)
    {
      pt1=pt2<-0
      for(v in 1:nV[i])
      {
        u1<-(t(X[[i]][v,])%*%c(beta0.vec.new[i],beta1.new))[1,1]
        pt1<-pt1+(1-score[i,v]*exp(-u1))*nT[i,v]/2
        pt2<-pt2+(score[i,v]*exp(-u1))*nT[i,v]/2
      }
      pt1<-pt1+(beta0.vec.new[i]-beta0.new)/sigma2.new
      pt2<-pt2+1/sigma2.new
      
      beta0.vec.upd[i]<-beta0.vec.new[i]-(pt1/pt2)
    }
    # update fix beta0
    beta0.upd<-mean(beta0.vec.upd)
    # update beta0 variance
    sigma2.upd<-mean((beta0.vec.upd-beta0.upd)^2)
    
    # update beta1
    pt1<-matrix(0,nrow=q-1,ncol=1)
    pt2<-matrix(0,nrow=q-1,ncol=q-1)
    for(i in 1:n)
    {
      for(v in 1:nV[i])
      {
        u1<-(t(X[[i]][v,])%*%c(beta0.vec.upd[i],beta1.new))[1,1]
        pt1<-pt1+(as.matrix(X[[i]][v,-1])-score[i,v]*exp(-u1)*as.matrix(X[[i]][v,-1]))*nT[i,v]/2
        pt2<-pt2+(score[i,v]*exp(-u1)*as.matrix(X[[i]][v,-1])%*%t(as.matrix(X[[i]][v,-1])))*nT[i,v]/2
      }
    }
    beta1.upd<-beta1.new-ginv(pt2)%*%pt1
    
    # calculate converge criterion
    diff<-max(abs(c(beta0.upd-beta0.new,beta1.upd-beta1.new)))
    
    # update parameters
    beta0.vec.new<-beta0.vec.upd
    beta0.new<-beta0.upd
    sigma2.new<-sigma2.upd
    beta1.new<-beta1.upd
    
    if(trace)
    {
      beta1.trace<-cbind(beta1.trace,beta1.new)
      beta0.vec.trace<-cbind(beta0.vec.trace,beta0.vec.new)
      beta0.trace<-cbind(beta0.trace,beta0.new)
      sigma2.trace<-cbind(sigma2.trace,sigma2.new)
      
      obj<-c(obj,obj.func(X,nT,S.iv,gamma,beta0.vec.new,beta0.new,sigma2.new,beta1.new))
    }
    
    # print(c(diff,obj.func(X,nT,S.iv,gamma,beta0.vec.new,beta0.new,sigma2.new,beta1.new)))
  }
  
  rownames(score)<-names(X)
  colnames(score)<-paste0("visit",1:max(nV))
  
  if(trace)
  {
    colnames(beta1.trace)=colnames(beta0.vec.trace)=colnames(beta0.trace)=colnames(sigma2.trace)<-paste0("iteration",0:(ncol(beta1.trace)-1))
    
    beta.trace<-rbind(beta0.trace,beta1.trace)
    rownames(beta.trace)<-colnames(X[[1]])
    
    if(score.return)
    {
      re<-list(gamma=c(gamma),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),score=score,
               beta.trace=beta.trace,beta0.random.trace=beta0.vec.trace,beta0.sigma2.trace=sigma2.trace,obj=obj)
    }else
    {
      re<-list(gamma=c(gamma),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
               beta.trace=beta.trace,beta0.random.trace=beta0.vec.trace,beta0.sigma2.trace=sigma2.trace,obj=obj)
    }
  }else
  {
    if(score.return)
    {
      re<-list(gamma=c(gamma),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),score=score)
    }else
    {
      re<-list(gamma=c(gamma),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr))
    }
  }
  
  return(re)
}
cap_beta<-function(Y,X,gamma,cov.shrinkage=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE)
{
  # Y: outcome list (longitudinal)
  # X: covariates (longitudinal)
  # gamma: linear projection
  
  #--------------------------------------------
  n<-length(Y)
  nV<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])
  
  q<-ncol(X[[1]])
  
  nT<-matrix(NA,n,max(nV))
  for(i in 1:n)
  {
    for(v in 1:nV[i])
    {
      nT[i,v]<-nrow(Y[[i]][[v]])
    }
    
    if(is.null(colnames(X[[i]]))==TRUE)
    {
      colnames(X[[i]])<-paste0("X",0:(q-1))
    }
  }
  
  if(min(nT,na.rm=TRUE)-5<p)
  {
    cov.shrinkage<-TRUE
  }
  #--------------------------------------------
  
  #--------------------------------------------
  # Initial value: estimate covariance matrix for each subject
  S.iv<-vector("list",length=n)
  score<-matrix(NA,n,max(nV))
  for(i in 1:n)
  {
    S.iv[[i]]<-array(NA,c(p,p,nV[i]))
    for(v in 1:nV[i])
    {
      if(cov.shrinkage)
      {
        S.iv[[i]][,,v]<-cov.ls(Y[[i]][[v]])
      }else
      {
        # sample covariance matrix
        S.iv[[i]][,,v]<-cov(Y[[i]][[v]])*(nT[i,v]-1)/nT[i,v]
      }
      
      score[i,v]<-t(gamma)%*%S.iv[[i]][,,v]%*%gamma
    }
  }
  #--------------------------------------------
  
  #--------------------------------------------
  # inintial estimate using a mixed effects model
  dtmp<-NULL
  for(i in 1:n)
  {
    dtmp<-rbind(dtmp,data.frame(ID=rep(i,nV[i]),score=log(score[i,1:nV[i]]),X[[i]][,-1]))
  }
  colnames(dtmp)<-c("ID","score",colnames(X[[1]])[-1])
  eval(parse(text=paste0("fit.tmp<-lme(score~",paste(colnames(X[[1]])[-1],collapse="+"),",random=~1|ID,data=dtmp,control=lmeControl(opt='optim'))")))
  coef.fix<-fit.tmp$coefficients$fixed
  beta1.new<-as.matrix(coef.fix[-1])
  beta0.vec.new<-coef.fix[1]+fit.tmp$coefficients$random$ID
  beta0.new<-mean(beta0.vec.new)
  sigma2.new<-mean((beta0.vec.new-beta0.new)^2)
  #--------------------------------------------
  
  if(trace)
  {
    beta1.trace<-cbind(beta1.new)
    beta0.vec.trace<-cbind(beta0.vec.new)
    beta0.trace<-cbind(beta0.new)
    sigma2.trace<-cbind(sigma2.new)
    
    if(cov.shrinkage)
    {
      par.sk<-NULL
    }
    
    obj<-obj.func(X,nT,S.iv,gamma,beta0.vec.new,beta0.new,sigma2.new,beta1.new)
  }
  
  s<-0
  diff<-100
  while(s<=max.itr&diff>tol)
  {
    s<-s+1
    
    # update random beta0
    beta0.vec.upd<-rep(NA,n)
    for(i in 1:n)
    {
      pt1=pt2<-0
      for(v in 1:nV[i])
      {
        u1<-(t(X[[i]][v,])%*%c(beta0.vec.new[i],beta1.new))[1,1]
        pt1<-pt1+(1-score[i,v]*exp(-u1))*nT[i,v]/2
        pt2<-pt2+(score[i,v]*exp(-u1))*nT[i,v]/2
      }
      pt1<-pt1+(beta0.vec.new[i]-beta0.new)/sigma2.new
      pt2<-pt2+1/sigma2.new
      
      beta0.vec.upd[i]<-beta0.vec.new[i]-(pt1/pt2)
    }
    # update fix beta0
    beta0.upd<-mean(beta0.vec.upd)
    # update beta0 variance
    sigma2.upd<-mean((beta0.vec.upd-beta0.upd)^2)
    
    # update beta1
    pt1<-matrix(0,nrow=q-1,ncol=1)
    pt2<-matrix(0,nrow=q-1,ncol=q-1)
    for(i in 1:n)
    {
      for(v in 1:nV[i])
      {
        u1<-(t(X[[i]][v,])%*%c(beta0.vec.upd[i],beta1.new))[1,1]
        pt1<-pt1+(as.matrix(X[[i]][v,-1])-score[i,v]*exp(-u1)*as.matrix(X[[i]][v,-1]))*nT[i,v]/2
        pt2<-pt2+(score[i,v]*exp(-u1)*as.matrix(X[[i]][v,-1])%*%t(as.matrix(X[[i]][v,-1])))*nT[i,v]/2
      }
    }
    beta1.upd<-beta1.new-ginv(pt2)%*%pt1
    
    if(cov.shrinkage)
    {
      # update shrinkage estimator of covariance matrix
      cov.upd<-cov.ls.const(Y,X,gamma,beta0.vec.upd,beta1.upd)
      S.iv<-cov.upd$S
      score<-matrix(NA,n,max(nV))
      for(i in 1:n)
      {
        for(v in 1:nV[i])
        {
          score[i,v]<-t(gamma)%*%S.iv[[i]][,,v]%*%gamma
        }
      }
    }
    
    # calculate converge criterion
    diff<-max(abs(c(beta0.upd-beta0.new,beta1.upd-beta1.new)))
    
    # update parameters
    beta0.vec.new<-beta0.vec.upd
    beta0.new<-beta0.upd
    sigma2.new<-sigma2.upd
    beta1.new<-beta1.upd
    
    if(trace)
    {
      beta1.trace<-cbind(beta1.trace,beta1.new)
      beta0.vec.trace<-cbind(beta0.vec.trace,beta0.vec.new)
      beta0.trace<-cbind(beta0.trace,beta0.new)
      sigma2.trace<-cbind(sigma2.trace,sigma2.new)
      
      if(cov.shrinkage)
      {
        par.sk<-cbind(par.sk,c(cov.upd$mu,cov.upd$delta,cov.upd$psi,cov.upd$phi,cov.upd$rho1,cov.upd$rho2))
        rownames(par.sk)<-c("mu","delta","psi","phi","rho1","rho2")
      }
      
      obj<-c(obj,obj.func(X,nT,S.iv,gamma,beta0.vec.new,beta0.new,sigma2.new,beta1.new))
    }
    
    # print(c(diff,obj.func(X,nT,S.iv,gamma,beta0.vec.new,beta0.new,sigma2.new,beta1.new)))
  }
  
  rownames(score)<-names(X)
  colnames(score)<-paste0("visit",1:max(nV))
  
  if(trace)
  {
    colnames(beta1.trace)=colnames(beta0.vec.trace)=colnames(beta0.trace)=colnames(sigma2.trace)<-paste0("iteration",0:(ncol(beta1.trace)-1))
    if(cov.shrinkage)
    {
      colnames(par.sk)<-paste0("iteration",1:ncol(par.sk))
    }
    
    beta.trace<-rbind(beta0.trace,beta1.trace)
    rownames(beta.trace)<-colnames(X[[1]])
    
    if(score.return)
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=c(gamma),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
                 shrinkage=data.frame(rho1=cov.upd$rho1,rho2=cov.upd$rho2,mu=cov.upd$mu,phi2=cov.upd$phi,psi2=cov.upd$psi,delta2=cov.upd$delta),score=score,
                 beta.trace=beta.trace,beta0.random.trace=beta0.vec.trace,beta0.sigma2.trace=sigma2.trace,shrinkage.trace=par.sk,obj=obj)
      }else
      {
        re<-list(gamma=c(gamma),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),score=score,
                 beta.trace=beta.trace,beta0.random.trace=beta0.vec.trace,beta0.sigma2.trace=sigma2.trace,obj=obj)
      }
      
    }else
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=c(gamma),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
                 shrinkage=data.frame(rho1=cov.upd$rho1,rho2=cov.upd$rho2,mu=cov.upd$mu,phi2=cov.upd$phi,psi2=cov.upd$psi,delta2=cov.upd$delta),
                 beta.trace=beta.trace,beta0.random.trace=beta0.vec.trace,beta0.sigma2.trace=sigma2.trace,shrinkage.trace=par.sk,obj=obj)
      }else
      {
        re<-list(gamma=c(gamma),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
                 beta.trace=beta.trace,beta0.random.trace=beta0.vec.trace,beta0.sigma2.trace=sigma2.trace,obj=obj)
      }
    }
  }else
  {
    if(score.return)
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=c(gamma),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
                 shrinkage=data.frame(rho1=cov.upd$rho1,rho2=cov.upd$rho2,mu=cov.upd$mu,phi2=cov.upd$phi,psi2=cov.upd$psi,delta2=cov.upd$delta),score=score)
      }else
      {
        re<-list(gamma=c(gamma),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),score=score)
      }
      
    }else
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=c(gamma),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
                 shrinkage=data.frame(rho1=cov.upd$rho1,rho2=cov.upd$rho2,mu=cov.upd$mu,phi2=cov.upd$phi,psi2=cov.upd$psi,delta2=cov.upd$delta))
      }else
      {
        re<-list(gamma=c(gamma),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr))
      }
    }
  }
  
  return(re)
}
#################################################

#################################################
# estimate both gamma and beta
# finding first direction
gamma.solve<-function(A,H)
{
  p<-ncol(H)
  
  H.svd<-svd(H)
  H.d.sqrt<-diag(sqrt(H.svd$d))
  H.d.sqrt.inv<-diag(1/sqrt(H.svd$d))
  H.sqrt.inv<-H.svd$u%*%H.d.sqrt.inv%*%t(H.svd$v)
  
  #---------------------------------------------------
  # svd decomposition method
  
  eigen.tmp<-eigen(H.d.sqrt.inv%*%t(H.svd$u)%*%A%*%H.svd$u%*%H.d.sqrt.inv)
  eigen.tmp.vec<-Re(eigen.tmp$vectors)
  re<-H.svd$u%*%H.d.sqrt.inv%*%eigen.tmp.vec[,p]
  
  # obj<-rep(NA,ncol(eigen.tmp$vectors))
  # for(j in 1:ncol(eigen.tmp$vectors))
  # {
  #   otmp<-H.svd$u%*%H.d.sqrt.inv%*%eigen.tmp.vec[,j]
  #   obj[j]<-t(otmp)%*%A%*%otmp
  # }
  # re<-H.svd$u%*%H.d.sqrt.inv%*%eigen.tmp.vec[,which.min(obj)]
  #---------------------------------------------------
  
  #---------------------------------------------------
  # eigenvector of A with respect to H
  
  # xtmp<-H.sqrt.inv%*%Re(eigen(H.sqrt.inv%*%A%*%H.sqrt.inv)$vectors)
  # opt.idx<-which.min(diag(t(xtmp)%*%A%*%xtmp))
  # re<-xtmp[,opt.idx]
  # re<-xtmp[,opt.idx]/sqrt(sum((xtmp[,opt.idx])^2))
  
  # eigen.tmp<-eigen(H.sqrt.inv%*%A%*%H.sqrt.inv)
  # eigen.tmp.vec<-Re(eigen.tmp$vectors)
  # 
  # obj<-rep(NA,ncol(eigen.tmp$vectors))
  # for(j in 1:ncol(eigen.tmp$vectors))
  # {
  #   otmp<-H.sqrt.inv%*%eigen.tmp.vec[,j]
  #   obj[j]<-t(otmp)%*%A%*%otmp
  # }
  # re<-H.sqrt.inv%*%eigen.tmp.vec[,p]
  #---------------------------------------------------
  
  return(re)
}

# estimate both gamma and beta
cap.cov_D1<-function(Y.cov,X,Y.n,method=c("CAP"),max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0=NULL)
{
  # Y.cov: a list of covariance matrix of Y
  # X: covariates (longitudinal)
  # Y.n: n by nV matrix of sample size
  
  #--------------------------------------------
  n<-length(Y.cov)
  nV<-sapply(X,nrow)
  p<-dim(Y.cov[[1]])[1]
  
  q<-ncol(X[[1]])
  
  nT<-Y.n
  for(i in 1:n)
  {
    if(is.null(colnames(X[[i]]))==TRUE)
    {
      colnames(X[[i]])<-paste0("X",0:(q-1))
    }
  }
  #--------------------------------------------
  
  #--------------------------------------------
  # Initial value: estimate covariance matrix for each subject
  S.iv<-Y.cov
  H<-matrix(0,p,p)
  for(i in 1:n)
  {
    for(v in 1:nV[i])
    {
      # H: average sample covariance matrix
      H<-H+S.iv[[i]][,,v]*(nT[i,v]/sum(nT,na.rm=TRUE))
    }
  }
  #--------------------------------------------
  
  #============================================
  if(method[1]=="CAP")
  {
    #--------------------------------------------
    # set initial value of gamma0
    if(is.null(gamma0))
    {
      set.seed(100)
      gamma.tmp<-matrix(rnorm((p+1+5)*p,mean=0,sd=1),nrow=p)
      gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      gamma0<-gamma0.mat[,sample(ncol(gamma0.mat),1)]
    }
    
    # calculate score
    score<-matrix(NA,n,max(nV))
    for(i in 1:n)
    {
      for(v in 1:nV[i])
      {
        score[i,v]<-t(gamma0)%*%S.iv[[i]][,,v]%*%gamma0
      }
    }
    
    # inintial estimate of beta's using a mixed effects model
    dtmp<-NULL
    for(i in 1:n)
    {
      dtmp<-rbind(dtmp,data.frame(ID=rep(i,nV[i]),score=log(score[i,1:nV[i]]),X[[i]][,-1]))
    }
    colnames(dtmp)<-c("ID","score",colnames(X[[1]])[-1])
    eval(parse(text=paste0("fit.tmp<-lme(score~",paste(colnames(X[[1]])[-1],collapse="+"),",random=~1|ID,data=dtmp,control=lmeControl(opt='optim'))")))
    coef.fix<-fit.tmp$coefficients$fixed
    beta1.new<-as.matrix(coef.fix[-1])
    beta0.vec.new<-coef.fix[1]+fit.tmp$coefficients$random$ID
    beta0.new<-mean(beta0.vec.new)
    sigma2.new<-mean((beta0.vec.new-beta0.new)^2)
    #--------------------------------------------
    
    if(trace)
    {
      gamma.trace<-cbind(gamma0)
      
      beta1.trace<-cbind(beta1.new)
      beta0.vec.trace<-cbind(beta0.vec.new)
      beta0.trace<-cbind(beta0.new)
      sigma2.trace<-cbind(sigma2.new)
      
      obj<-obj.func(X,nT,S.iv,gamma0,beta0.vec.new,beta0.new,sigma2.new,beta1.new)
    }
    
    s<-0
    diff<-100
    while(s<=max.itr&diff>tol)
    {
      s<-s+1
      
      score<-matrix(NA,n,max(nV))
      for(i in 1:n)
      {
        for(v in 1:nV[i])
        {
          score[i,v]<-t(gamma0)%*%S.iv[[i]][,,v]%*%gamma0
        }
      }
      
      #--------------------------------------------
      # update random beta0
      beta0.vec.upd<-rep(NA,n)
      for(i in 1:n)
      {
        pt1=pt2<-0
        for(v in 1:nV[i])
        {
          u1<-(t(X[[i]][v,])%*%c(beta0.vec.new[i],beta1.new))[1,1]
          pt1<-pt1+(1-score[i,v]*exp(-u1))*nT[i,v]/2
          pt2<-pt2+(score[i,v]*exp(-u1))*nT[i,v]/2
        }
        pt1<-pt1+(beta0.vec.new[i]-beta0.new)/sigma2.new
        pt2<-pt2+1/sigma2.new
        
        beta0.vec.upd[i]<-beta0.vec.new[i]-(pt1/pt2)
      }
      # update fix beta0
      beta0.upd<-mean(beta0.vec.upd)
      # update beta0 variance
      sigma2.upd<-mean((beta0.vec.upd-beta0.upd)^2)
      
      # update beta1
      pt1<-matrix(0,nrow=q-1,ncol=1)
      pt2<-matrix(0,nrow=q-1,ncol=q-1)
      for(i in 1:n)
      {
        for(v in 1:nV[i])
        {
          u1<-(t(X[[i]][v,])%*%c(beta0.vec.upd[i],beta1.new))[1,1]
          pt1<-pt1+(as.matrix(X[[i]][v,-1])-score[i,v]*exp(-u1)*as.matrix(X[[i]][v,-1]))*nT[i,v]/2
          pt2<-pt2+(score[i,v]*exp(-u1)*as.matrix(X[[i]][v,-1])%*%t(as.matrix(X[[i]][v,-1])))*nT[i,v]/2
        }
      }
      beta1.upd<-beta1.new-ginv(pt2)%*%pt1
      
      # update gamma
      Amat<-matrix(0,p,p)
      for(i in 1:n)
      {
        for(v in 1:nV[i])
        {
          u1<-(t(X[[i]][v,])%*%c(beta0.vec.new[i],beta1.new))[1,1]
          Amat<-Amat+exp(-u1)*S.iv[[i]][,,v]*nT[i,v]/2
        }
      }
      gamma.upd<-gamma.solve(Amat,H)
      #--------------------------------------------
      
      # calculate converge criterion
      # diff<-max(c(max(abs(gamma.upd-gamma0)),max(abs(beta0.upd-beta0.new)),max(abs(beta1.upd-beta1.new))))
      diff<-max(abs(c(beta0.upd-beta0.new,beta1.upd-beta1.new)))
      
      # update parameters
      beta0.vec.new<-beta0.vec.upd
      beta0.new<-beta0.upd
      sigma2.new<-sigma2.upd
      beta1.new<-beta1.upd
      gamma0<-gamma.upd
      
      if(trace)
      {
        gamma.trace<-cbind(gamma.trace,gamma0)
        
        beta1.trace<-cbind(beta1.trace,beta1.new)
        beta0.vec.trace<-cbind(beta0.vec.trace,beta0.vec.new)
        beta0.trace<-cbind(beta0.trace,beta0.new)
        sigma2.trace<-cbind(sigma2.trace,sigma2.new)

        obj<-c(obj,obj.func(X,nT,S.iv,gamma0,beta0.vec.new,beta0.new,sigma2.new,beta1.new))
      }
      
      # print(c(diff,obj.func(X,nT,S.iv,gamma0,beta0.vec.new,beta0.new,sigma2.new,beta1.new)))
    }
  }
  #============================================
  
  # standardize gamma estimate
  gamma0<-c(gamma0)/sqrt(sum(gamma0^2))
  if(gamma0[which.max(abs(gamma0))]<0)
  {
    gamma0<--gamma0
  }
  # reestimate beta
  beta.out<-cap.cov_beta(Y.cov,X,Y.n,gamma0,max.itr=max.itr,tol=tol,score.return=TRUE,trace=FALSE)
  beta0.vec.new<-beta.out$beta0.random
  beta0.new<-beta.out$beta[1]
  beta1.new<-beta.out$beta[-1]
  sigma2.new<-beta.out$beta0.sigma2
  
  score<-beta.out$score
  rownames(score)<-names(X)
  colnames(score)<-paste0("visit",1:max(nV))
  
  if(trace)
  {
    colnames(beta1.trace)=colnames(beta0.vec.trace)=colnames(beta0.trace)=colnames(sigma2.trace)=colnames(gamma.trace)<-paste0("iteration",0:(ncol(beta1.trace)-1))
    
    beta.trace<-rbind(beta0.trace,beta1.trace)
    rownames(beta.trace)<-colnames(X[[1]])
    
    if(score.return)
    {
      re<-list(gamma=c(gamma0),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),score=score,
               beta.trace=beta.trace,beta0.random.trace=beta0.vec.trace,beta0.sigma2.trace=sigma2.trace,obj=obj)
    }else
    {
      re<-list(gamma=c(gamma0),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
               beta.trace=beta.trace,beta0.random.trace=beta0.vec.trace,beta0.sigma2.trace=sigma2.trace,obj=obj)
    }
  }else
  {
    if(score.return)
    {
      re<-list(gamma=c(gamma0),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),score=score)
    }else
    {
      re<-list(gamma=c(gamma0),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr))
    }
  }
  
  return(re)
}
cap_D1<-function(Y,X,method=c("CAP"),cov.shrinkage=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0=NULL)
{
  # Y: outcome list (longitudinal)
  # X: covariates (longitudinal)
  
  #--------------------------------------------
  n<-length(Y)
  nV<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])
  
  q<-ncol(X[[1]])
  
  nT<-matrix(NA,n,max(nV))
  for(i in 1:n)
  {
    for(v in 1:nV[i])
    {
      nT[i,v]<-nrow(Y[[i]][[v]])
    }
    
    if(is.null(colnames(X[[i]]))==TRUE)
    {
      colnames(X[[i]])<-paste0("X",0:(q-1))
    }
  }
  
  if(min(nT,na.rm=TRUE)-5<p)
  {
    cov.shrinkage<-TRUE
  }
  #--------------------------------------------
  
  #--------------------------------------------
  # Initial value: estimate covariance matrix for each subject
  S.iv<-vector("list",length=n)
  H<-matrix(0,p,p)
  for(i in 1:n)
  {
    S.iv[[i]]<-array(NA,c(p,p,nV[i]))
    for(v in 1:nV[i])
    {
      if(cov.shrinkage)
      {
        S.iv[[i]][,,v]<-cov.ls(Y[[i]][[v]])
      }else
      {
        # sample covariance matrix
        S.iv[[i]][,,v]<-cov(Y[[i]][[v]])*(nT[i,v]-1)/nT[i,v]
      }
      # H: average sample covariance matrix
      H<-H+S.iv[[i]][,,v]*(nT[i,v]/sum(nT,na.rm=TRUE))
    }
  }
  #--------------------------------------------
  
  #============================================
  if(method[1]=="CAP")
  {
    #--------------------------------------------
    # set initial value of gamma0
    if(is.null(gamma0))
    {
      set.seed(100)
      gamma.tmp<-matrix(rnorm((p+1+5)*p,mean=0,sd=1),nrow=p)
      gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      gamma0<-gamma0.mat[,sample(ncol(gamma0.mat),1)]
    }
    
    # calculate score
    score<-matrix(NA,n,max(nV))
    for(i in 1:n)
    {
      for(v in 1:nV[i])
      {
        score[i,v]<-t(gamma0)%*%S.iv[[i]][,,v]%*%gamma0
      }
    }
    
    # inintial estimate of beta's using a mixed effects model
    dtmp<-NULL
    for(i in 1:n)
    {
      dtmp<-rbind(dtmp,data.frame(ID=rep(i,nV[i]),score=log(score[i,1:nV[i]]),X[[i]][,-1]))
    }
    colnames(dtmp)<-c("ID","score",colnames(X[[1]])[-1])
    eval(parse(text=paste0("fit.tmp<-lme(score~",paste(colnames(X[[1]])[-1],collapse="+"),",random=~1|ID,data=dtmp,control=lmeControl(opt='optim'))")))
    coef.fix<-fit.tmp$coefficients$fixed
    beta1.new<-as.matrix(coef.fix[-1])
    beta0.vec.new<-coef.fix[1]+fit.tmp$coefficients$random$ID
    beta0.new<-mean(beta0.vec.new)
    sigma2.new<-mean((beta0.vec.new-beta0.new)^2)
    #--------------------------------------------
    
    if(trace)
    {
      gamma.trace<-cbind(gamma0)
      
      beta1.trace<-cbind(beta1.new)
      beta0.vec.trace<-cbind(beta0.vec.new)
      beta0.trace<-cbind(beta0.new)
      sigma2.trace<-cbind(sigma2.new)
      
      if(cov.shrinkage)
      {
        par.sk<-NULL
      }
      
      obj<-obj.func(X,nT,S.iv,gamma0,beta0.vec.new,beta0.new,sigma2.new,beta1.new)
    }
    
    s<-0
    diff<-100
    while(s<=max.itr&diff>tol)
    {
      s<-s+1
      
      # calculate score
      score<-matrix(NA,n,max(nV))
      for(i in 1:n)
      {
        for(v in 1:nV[i])
        {
          score[i,v]<-t(gamma0)%*%S.iv[[i]][,,v]%*%gamma0
        }
      }
      
      #--------------------------------------------
      # update random beta0
      beta0.vec.upd<-rep(NA,n)
      for(i in 1:n)
      {
        pt1=pt2<-0
        for(v in 1:nV[i])
        {
          u1<-(t(X[[i]][v,])%*%c(beta0.vec.new[i],beta1.new))[1,1]
          pt1<-pt1+(1-score[i,v]*exp(-u1))*nT[i,v]/2
          pt2<-pt2+(score[i,v]*exp(-u1))*nT[i,v]/2
        }
        pt1<-pt1+(beta0.vec.new[i]-beta0.new)/sigma2.new
        pt2<-pt2+1/sigma2.new
        
        beta0.vec.upd[i]<-beta0.vec.new[i]-(pt1/pt2)
      }
      # update fix beta0
      beta0.upd<-mean(beta0.vec.upd)
      # update beta0 variance
      sigma2.upd<-mean((beta0.vec.upd-beta0.upd)^2)
      
      # update beta1
      pt1<-matrix(0,nrow=q-1,ncol=1)
      pt2<-matrix(0,nrow=q-1,ncol=q-1)
      for(i in 1:n)
      {
        for(v in 1:nV[i])
        {
          u1<-(t(X[[i]][v,])%*%c(beta0.vec.upd[i],beta1.new))[1,1]
          pt1<-pt1+(as.matrix(X[[i]][v,-1])-score[i,v]*exp(-u1)*as.matrix(X[[i]][v,-1]))*nT[i,v]/2
          pt2<-pt2+(score[i,v]*exp(-u1)*as.matrix(X[[i]][v,-1])%*%t(as.matrix(X[[i]][v,-1])))*nT[i,v]/2
        }
      }
      beta1.upd<-beta1.new-ginv(pt2)%*%pt1

      if(cov.shrinkage)
      {
        # update shrinkage estimator of covariance matrix
        cov.upd<-cov.ls.const(Y,X,gamma0,beta0.vec.upd,beta1.upd)
        S.iv<-cov.upd$S
        score<-matrix(NA,n,max(nV))
        H<-matrix(0,p,p)
        for(i in 1:n)
        {
          for(v in 1:nV[i])
          {
            score[i,v]<-t(gamma0)%*%S.iv[[i]][,,v]%*%gamma0
            
            # H: average sample covariance matrix
            H<-H+S.iv[[i]][,,v]*(nT[i,v]/sum(nT,na.rm=TRUE))
          }
        }
      }
      
      # update gamma
      Amat<-matrix(0,p,p)
      for(i in 1:n)
      {
        for(v in 1:nV[i])
        {
          u1<-(t(X[[i]][v,])%*%c(beta0.vec.new[i],beta1.new))[1,1]
          Amat<-Amat+exp(-u1)*S.iv[[i]][,,v]*nT[i,v]/2
        }
      }
      gamma.upd<-gamma.solve(Amat,H)
      #--------------------------------------------
      
      # calculate converge criterion
      # diff<-max(c(max(abs(gamma.upd-gamma0)),max(abs(beta0.upd-beta0.new)),max(abs(beta1.upd-beta1.new))))
      diff<-max(abs(c(beta0.upd-beta0.new,beta1.upd-beta1.new)))
      
      # update parameters
      beta0.vec.new<-beta0.vec.upd
      beta0.new<-beta0.upd
      sigma2.new<-sigma2.upd
      beta1.new<-beta1.upd
      gamma0<-gamma.upd
      
      if(trace)
      {
        gamma.trace<-cbind(gamma.trace,gamma0)
        
        beta1.trace<-cbind(beta1.trace,beta1.new)
        beta0.vec.trace<-cbind(beta0.vec.trace,beta0.vec.new)
        beta0.trace<-cbind(beta0.trace,beta0.new)
        sigma2.trace<-cbind(sigma2.trace,sigma2.new)
        
        if(cov.shrinkage)
        {
          par.sk<-cbind(par.sk,c(cov.upd$mu,cov.upd$delta,cov.upd$psi,cov.upd$phi,cov.upd$rho1,cov.upd$rho2))
          rownames(par.sk)<-c("mu","delta","psi","phi","rho1","rho2")
        }
        
        obj<-c(obj,obj.func(X,nT,S.iv,gamma0,beta0.vec.new,beta0.new,sigma2.new,beta1.new))
      }
      
      # print(c(diff,obj.func(X,nT,S.iv,gamma0,beta0.vec.new,beta0.new,sigma2.new,beta1.new)))
    }
  }
  #============================================
  
  # standardize gamma estimate
  gamma0<-c(gamma0)/sqrt(sum(gamma0^2))
  if(gamma0[which.max(abs(gamma0))]<0)
  {
    gamma0<--gamma0
  }
  # reestimate beta
  beta.out<-cap_beta(Y,X,gamma0,cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,score.return=TRUE,trace=FALSE)
  beta0.vec.new<-beta.out$beta0.random
  beta0.new<-beta.out$beta[1]
  beta1.new<-beta.out$beta[-1]
  sigma2.new<-beta.out$beta0.sigma2
  
  score<-beta.out$score
  rownames(score)<-names(X)
  colnames(score)<-paste0("visit",1:max(nV))
  
  if(cov.shrinkage)
  {
    cov.upd<-cov.ls.const(Y,X,gamma0,beta0.vec.upd,beta1.upd)
  }
  
  if(trace)
  {
    colnames(beta1.trace)=colnames(beta0.vec.trace)=colnames(beta0.trace)=colnames(sigma2.trace)=colnames(gamma.trace)<-paste0("iteration",0:(ncol(beta1.trace)-1))
    if(cov.shrinkage)
    {
      colnames(par.sk)<-paste0("iteration",1:ncol(par.sk))
    }
    
    beta.trace<-rbind(beta0.trace,beta1.trace)
    rownames(beta.trace)<-colnames(X[[1]])
    
    if(score.return)
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=c(gamma0),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
                 shrinkage=data.frame(rho1=cov.upd$rho1,rho2=cov.upd$rho2,mu=cov.upd$mu,phi2=cov.upd$phi,psi2=cov.upd$psi,delta2=cov.upd$delta),score=score,
                 beta.trace=beta.trace,beta0.random.trace=beta0.vec.trace,beta0.sigma2.trace=sigma2.trace,shrinkage.trace=par.sk,obj=obj)
      }else
      {
        re<-list(gamma=c(gamma0),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),score=score,
                 beta.trace=beta.trace,beta0.random.trace=beta0.vec.trace,beta0.sigma2.trace=sigma2.trace,obj=obj)
      }
      
    }else
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=c(gamma0),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
                 shrinkage=data.frame(rho1=cov.upd$rho1,rho2=cov.upd$rho2,mu=cov.upd$mu,phi2=cov.upd$phi,psi2=cov.upd$psi,delta2=cov.upd$delta),
                 beta.trace=beta.trace,beta0.random.trace=beta0.vec.trace,beta0.sigma2.trace=sigma2.trace,shrinkage.trace=par.sk,obj=obj)
      }else
      {
        re<-list(gamma=c(gamma0),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
                 beta.trace=beta.trace,beta0.random.trace=beta0.vec.trace,beta0.sigma2.trace=sigma2.trace,obj=obj)
      }
    }
  }else
  {
    if(score.return)
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=c(gamma0),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
                 shrinkage=data.frame(rho1=cov.upd$rho1,rho2=cov.upd$rho2,mu=cov.upd$mu,phi2=cov.upd$phi,psi2=cov.upd$psi,delta2=cov.upd$delta),score=score)
      }else
      {
        re<-list(gamma=c(gamma0),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),score=score)
      }
      
    }else
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=c(gamma0),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr),
                 shrinkage=data.frame(rho1=cov.upd$rho1,rho2=cov.upd$rho2,mu=cov.upd$mu,phi2=cov.upd$phi,psi2=cov.upd$psi,delta2=cov.upd$delta))
      }else
      {
        re<-list(gamma=c(gamma0),beta=c(beta0.new,beta1.new),beta0.random=beta0.vec.new,beta0.sigma2=sigma2.new,convergence=(s<max.itr))
      }
    }
  }
  
  return(re)
}

# try several initial value of gamma and optimize over the objective function
cap.cov_D1_opt<-function(Y.cov,X,Y.n,method=c("CAP"),max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,seed=500)
{
  # Y.cov: a list of covariance matrix of Y
  # X: covariates (longitudinal)
  # Y.n: n by nV matrix of sample size
  
  #--------------------------------------------
  n<-length(Y.cov)
  nV<-sapply(X,nrow)
  p<-dim(Y.cov[[1]])[1]
  
  q<-ncol(X[[1]])
  
  nT<-Y.n
  for(i in 1:n)
  {
    if(is.null(colnames(X[[i]]))==TRUE)
    {
      colnames(X[[i]])<-paste0("X",0:(q-1))
    }
  }
  #--------------------------------------------
  
  #--------------------------------------------
  # Initial value: estimate covariance matrix for each subject
  S.iv<-Y.cov
  H<-matrix(0,p,p)
  for(i in 1:n)
  {
    for(v in 1:nV[i])
    {
      # H: average sample covariance matrix
      H<-H+S.iv[[i]][,,v]*(nT[i,v]/sum(nT,na.rm=TRUE))
    }
  }
  #--------------------------------------------
  
  #============================================
  if(method[1]=="CAP")
  {
    # set initial values
    if(is.null(gamma0.mat))
    {
      #--------------------------------
      # gamma0.mat<-matrix(NA,p,p+1+5)
      # for(j in 1:p)
      # {
      #   gamma0.mat[,j]<-rep(0,p)
      #   gamma0.mat[j,j]<-1
      # }
      # gamma0.mat[,p+1]<-rep(1,p)/sqrt(sum(rep(1,p)^2))
      # 
      # set.seed(500)
      # gamma.tmp<-matrix(rnorm(5*p,mean=0,sd=1),nrow=p)
      # gamma0.mat[,(p+2):(p+1+5)]<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      #--------------------------------
      
      #--------------------------------
      set.seed(seed)
      gamma.tmp<-matrix(rnorm((p+1+5)*p,mean=0,sd=1),nrow=p)
      gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      #--------------------------------
    }
    if(is.null(ninitial))
    {
      ninitial<-min(ncol(gamma0.mat),10)
    }else
    {
      if(ninitial>ncol(gamma0.mat))
      {
        ninitial<-ncol(gamma0.mat)
      }
    }
    set.seed(seed)
    gamma0.mat<-matrix(gamma0.mat[,sort(sample(1:ncol(gamma0.mat),ninitial,replace=FALSE))],ncol=ninitial)
    
    re.tmp<-vector("list",ncol(gamma0.mat))
    obj<-rep(NA,ncol(gamma0.mat))
    for(kk in 1:ncol(gamma0.mat))
    {
      try(re.tmp[[kk]]<-cap.cov_D1(Y.cov,X,Y.n,method=method[1],max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0=gamma0.mat[,kk]))
      
      if(is.null(re.tmp[[kk]])==FALSE)
      {
        #--------------------------------
        gamma.unscale<-re.tmp[[kk]]$gamma/sqrt(t(re.tmp[[kk]]$gamma)%*%H%*%re.tmp[[kk]]$gamma)[1,1]
        try(beta.tmp<-cap.cov_beta(Y.cov,X,Y.n,gamma.unscale,max.itr=max.itr,tol=tol,trace=FALSE,score.return=FALSE))
        try(obj[kk]<-obj.func(X,nT,Y.cov,gamma.unscale,beta.tmp$beta0.random,beta.tmp$beta[1],beta.tmp$beta0.sigma2,beta.tmp$beta[-1]))
        #--------------------------------
      }
    }
    opt.idx<-which.min(obj)
    re<-re.tmp[[opt.idx]]
  }
  #============================================
  
  return(re)
}
cap_D1_opt<-function(Y,X,method=c("CAP"),cov.shrinkage=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,seed=500)
{
  # Y: outcome list (longitudinal)
  # X: covariates (longitudinal)
  
  #--------------------------------------------
  n<-length(Y)
  nV<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])
  
  q<-ncol(X[[1]])
  
  nT<-matrix(NA,n,max(nV))
  for(i in 1:n)
  {
    for(v in 1:nV[i])
    {
      nT[i,v]<-nrow(Y[[i]][[v]])
    }
    
    if(is.null(colnames(X[[i]]))==TRUE)
    {
      colnames(X[[i]])<-paste0("X",0:(q-1))
    }
  }
  
  if(min(nT,na.rm=TRUE)-5<p)
  {
    cov.shrinkage<-TRUE
  }
  #--------------------------------------------
  
  #--------------------------------------------
  # Initial value: estimate covariance matrix for each subject
  S.iv<-vector("list",length=n)
  H<-matrix(0,p,p)
  for(i in 1:n)
  {
    S.iv[[i]]<-array(NA,c(p,p,nV[i]))
    for(v in 1:nV[i])
    {
      if(cov.shrinkage)
      {
        S.iv[[i]][,,v]<-cov.ls(Y[[i]][[v]])
      }else
      {
        # sample covariance matrix
        S.iv[[i]][,,v]<-cov(Y[[i]][[v]])*(nT[i,v]-1)/nT[i,v]
      }
      # H: average sample covariance matrix
      H<-H+S.iv[[i]][,,v]*(nT[i,v]/sum(nT,na.rm=TRUE))
    }
  }
  #--------------------------------------------
  
  #============================================
  if(method[1]=="CAP")
  {
    # set initial values
    if(is.null(gamma0.mat))
    {
      #--------------------------------
      # gamma0.mat<-matrix(NA,p,p+1+5)
      # for(j in 1:p)
      # {
      #   gamma0.mat[,j]<-rep(0,p)
      #   gamma0.mat[j,j]<-1
      # }
      # gamma0.mat[,p+1]<-rep(1,p)/sqrt(sum(rep(1,p)^2))
      # 
      # set.seed(500)
      # gamma.tmp<-matrix(rnorm(5*p,mean=0,sd=1),nrow=p)
      # gamma0.mat[,(p+2):(p+1+5)]<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      #--------------------------------
      
      #--------------------------------
      set.seed(seed)
      gamma.tmp<-matrix(rnorm((p+1+5)*p,mean=0,sd=1),nrow=p)
      gamma0.mat<-apply(gamma.tmp,2,function(x){return(x/sqrt(sum(x^2)))})
      #--------------------------------
    }
    if(is.null(ninitial))
    {
      ninitial<-min(ncol(gamma0.mat),10)
    }else
    {
      if(ninitial>ncol(gamma0.mat))
      {
        ninitial<-ncol(gamma0.mat)
      }
    }
    set.seed(seed)
    gamma0.mat<-matrix(gamma0.mat[,sort(sample(1:ncol(gamma0.mat),ninitial,replace=FALSE))],ncol=ninitial)
    
    re.tmp<-vector("list",ncol(gamma0.mat))
    obj<-rep(NA,ncol(gamma0.mat))
    for(kk in 1:ncol(gamma0.mat))
    {
      try(re.tmp[[kk]]<-cap_D1(Y,X,method=method[1],cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0=gamma0.mat[,kk]))
      
      if(is.null(re.tmp[[kk]])==FALSE)
      {
        #--------------------------------
        gamma.unscale<-re.tmp[[kk]]$gamma/sqrt(t(re.tmp[[kk]]$gamma)%*%H%*%re.tmp[[kk]]$gamma)[1,1]
        try(beta.tmp<-cap_beta(Y,X,gamma.unscale,cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,trace=FALSE,score.return=FALSE))
        if(cov.shrinkage)
        {
          cov.est<-cov.ls.const(Y,X,gamma.unscale,beta.tmp$beta0.random,beta.tmp$beta[-1])
          try(obj[kk]<-obj.func(X,nT,cov.est$S,gamma.unscale,beta.tmp$beta0.random,beta.tmp$beta[1],beta.tmp$beta0.sigma2,beta.tmp$beta[-1]))
        }else
        {
          try(obj[kk]<-obj.func(X,nT,S.iv,gamma.unscale,beta.tmp$beta0.random,beta.tmp$beta[1],beta.tmp$beta0.sigma2,beta.tmp$beta[-1]))
        }
        #--------------------------------
      }
    }
    opt.idx<-which.min(obj)
    re<-re.tmp[[opt.idx]]
  }
  #============================================
  
  return(re)
}
#################################################

#################################################
# second and higher direction
cap_Dk<-function(Y,X,Phi0=NULL,method=c("CAP"),OC=FALSE,cov.shrinkage=TRUE,max.itr=1000,tol=1e-4,score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,seed=500)
{
  # Y: outcome list (longitudinal)
  # X: covariates (longitudinal)
  # Phi0: identified components
  
  if(is.null(Phi0))
  {
    return(cap_D1_opt(Y,X,method=method[1],cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed))
  }else
  {
    #--------------------------------------------
    n<-length(Y)
    nV<-sapply(Y,length)
    p<-ncol(Y[[1]][[1]])
    
    q<-ncol(X[[1]])
    
    nT<-matrix(NA,n,max(nV))
    for(i in 1:n)
    {
      for(v in 1:nV[i])
      {
        nT[i,v]<-nrow(Y[[i]][[v]])
      }
      
      if(is.null(colnames(X[[i]]))==TRUE)
      {
        colnames(X[[i]])<-paste0("X",0:(q-1))
      }
    }
    
    if(min(nT,na.rm=TRUE)-5<p)
    {
      cov.shrinkage<-TRUE
    }
    #--------------------------------------------
    
    p0<-ncol(Phi0)
    # estimate beta
    beta.est<-vector("list",length=p0)
    names(beta.est)<-paste0("D",1:p0)
    for(j in 1:p0)
    {
      beta.est[[j]]<-cap_beta(Y,X,gamma=Phi0[,j],cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,trace=FALSE)
    }
    Ytmp<-vector("list",length=n)
    for(i in 1:n)
    {
      Ytmp[[i]]<-vector("list",length=nV[i])
      for(v in 1:nV[i])
      {
        Y2tmp<-Y[[i]][[v]]-Y[[i]][[v]]%*%(Phi0%*%t(Phi0))
        if(cov.shrinkage)
        {
          #----------------------------------
          # no need to add back the intercept using shrinkage method
          Ytmp[[i]][[v]]<-Y2tmp
          #----------------------------------
        }else
        {
          beta0.tmp<-rep(NA,p0)
          for(j in 1:p0)
          {
            beta0.tmp[j]<-beta.est[[j]]$beta0.random[i]
          }
          #----------------------------------
          Y2tmp.svd<-svd(Y2tmp)
          Ytmp[[i]][[v]]<-Y2tmp.svd$u%*%diag(c(Y2tmp.svd$d[1:(p-p0)],sqrt(exp(beta0.tmp)*nT[i,v])))%*%t(Y2tmp.svd$v)
          #----------------------------------
        }
      }
    }
    
    if(method[1]=="CAP")
    {
      if(OC==FALSE)
      {
        re.tmp<-cap_D1_opt(Ytmp,X,method=method[1],cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE,gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed)
      }else
      {
        re.tmp<-cap_D1_opt(Ytmp,X,method=method[1],cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE,gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed) 
      }
    }
    
    re<-re.tmp
    
    re$orthogonal<-c(t(re.tmp$gamma)%*%Phi0)
    
    return(re)
  }
}
#################################################

#################################################
# level of diagonalization
diag.level<-function(Y,Phi,cov.shrinkage=TRUE)
{
  # Y: outcome list (longitudinal)
  # Phi: rotation matrix
  
  if(is.null(ncol(Phi))|ncol(Phi)==1)
  {
    return("Dimension of Phi is less than 2")
  }else
  {
    n<-length(Y)
    nV<-sapply(Y,length)
    p<-ncol(Y[[1]][[1]])
    
    nT<-matrix(NA,n,max(nV))
    for(i in 1:n)
    {
      for(v in 1:nV[i])
      {
        nT[i,v]<-nrow(Y[[i]][[v]])
      }
    }
    
    ps<-ncol(Phi)
    
    #--------------------------------------------
    # estimated covariance matrix
    S.iv<-vector("list",length=n)
    for(i in 1:n)
    {
      S.iv[[i]]<-array(NA,c(p,p,nV[i]))
      for(v in 1:nV[i])
      {
        if(cov.shrinkage)
        {
          S.iv[[i]][,,v]<-cov.ls(Y[[i]][[v]])
        }else
        {
          # sample covariance matrix
          S.iv[[i]][,,v]<-cov(Y[[i]][[v]])*(nT[i,v]-1)/nT[i,v]
        }
      }
    }
    #--------------------------------------------
    
    #--------------------------------------------
    dl.sub<-array(NA,c(n,max(nV),ps))
    dimnames(dl.sub)[[3]]<-paste0("Dim",1:ps)
    for(i in 1:n)
    {
      for(v in 1:nV[i])
      {
        dl.sub[i,v,1]<-1
        for(j in 2:ps)
        {
          phi.tmp<-Phi[,1:j]
          mat.tmp<-t(phi.tmp)%*%S.iv[[i]][,,v]%*%phi.tmp
          dl.sub[i,v,j]<-(det(diag(diag(mat.tmp)))/det(mat.tmp))^(nT[i,v]/sum(nT,na.rm=TRUE))
        }
      }
    }
    pmean<-apply(dl.sub,3,function(x){return(prod(x,na.rm=TRUE))})
    #--------------------------------------------
    
    re<-list(avg.level=pmean,sub.level=dl.sub)
    return(re)
  }
}
#################################################

#################################################
# CAP function: either specify the number of components or choose the number of components based on DfD
capReg<-function(Y,X,stop.crt=c("nD","DfD"),nD=NULL,DfD.thred=5,method=c("CAP"),OC=FALSE,cov.shrinkage=TRUE,max.itr=1000,tol=1e-4,
                 score.return=TRUE,trace=FALSE,gamma0.mat=NULL,ninitial=NULL,seed=500,verbose=TRUE)
{
  # Y: outcome list (longitudinal)
  # X: covariates (longitudinal)
  # stop.crt: stopping criterion, nD=# of directions, DfD=DfD threshold
  
  if(stop.crt[1]=="nD"&is.null(nD))
  {
    stop.crt<-"DfD"
  }
  
  #--------------------------------------------
  n<-length(Y)
  nV<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])
  
  q<-ncol(X[[1]])
  
  nT<-matrix(NA,n,max(nV))
  for(i in 1:n)
  {
    for(v in 1:nV[i])
    {
      nT[i,v]<-nrow(Y[[i]][[v]])
    }
    
    if(is.null(colnames(X[[i]]))==TRUE)
    {
      colnames(X[[i]])<-paste0("X",0:(q-1))
    }
  }
  #--------------------------------------------
  
  if(min(nT,na.rm=TRUE)-5<p)
  {
    cov.shrinkage<-TRUE
  }
  
  if(method[1]=="CAP")
  {
    #--------------------------------------------
    # First direction
    tm1<-system.time(re1<-cap_D1_opt(Y,X,method="CAP",cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,score.return=TRUE,trace=FALSE,gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed))
    
    Phi.est<-matrix(re1$gamma,ncol=1)
    beta.est<-matrix(re1$beta,ncol=1)
    beta0.random.est<-matrix(re1$beta0.random,ncol=1)
    beta0.sigma2.est<-matrix(re1$beta0.sigma2,ncol=1)
    rownames(beta0.sigma2.est)<-"beta0.sigma2"
    
    cp.time<-matrix(as.numeric(tm1[1:3]),ncol=1)
    rownames(cp.time)<-c("user","system","elapsed")
    
    if(verbose)
    {
      print(paste0("Component ",ncol(Phi.est)))
    }
    
    if(cov.shrinkage)
    {
      sk.out<-matrix(re1$shrinkage,nrow=1)
      colnames(sk.out)<-colnames(re1$shrinkage)
    }
    #--------------------------------------------
    
    if(stop.crt[1]=="nD")
    {
      if(score.return)
      {
        score<-array(NA,c(n,max(nV),nD))
        dimnames(score)[[2]]<-paste0("visit",1:max(nV))
        dimnames(score)[[3]]<-paste0("D",1:nD)
        
        score[,,1]<-re1$score
      }
      
      if(nD>1)
      {
        for(j in 2:nD)
        {
          re.tmp<-NULL
          try(tm.tmp<-system.time(re.tmp<-cap_Dk(Y,X,Phi0=Phi.est,method="CAP",OC=OC,cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE,
                                                 gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed)))
          if(is.null(re.tmp)==FALSE)
          {
            Phi.est<-cbind(Phi.est,re.tmp$gamma)
            beta.est<-cbind(beta.est,re.tmp$beta)
            beta0.random.est<-cbind(beta0.random.est,re.tmp$beta0.random)
            beta0.sigma2.est<-cbind(beta0.sigma2.est,re.tmp$beta0.sigma2)
            
            cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
            
            if(verbose)
            {
              print(paste0("Component ",ncol(Phi.est)))
            }
            
            if(cov.shrinkage)
            {
              sk.out<-rbind(sk.out,re.tmp$shrinkage)
            }
            if(score.return)
            {
              score[,,j]<-re.tmp$score
            }
          }else
          {
            break
          }
        }
        
        colnames(Phi.est)=colnames(beta.est)=colnames(beta0.random.est)=colnames(beta0.sigma2.est)<-paste0("D",1:ncol(Phi.est))
        rownames(Phi.est)<-paste0("V",1:p)
        rownames(beta.est)<-colnames(X[[1]])
        cp.time<-cbind(cp.time,apply(cp.time,1,sum))
        colnames(cp.time)<-c(paste0("D",1:ncol(Phi.est)),"Total")
        if(cov.shrinkage)
        {
          rownames(sk.out)<-paste0("D",1:ncol(Phi.est))
        }
        
        if(ncol(Phi.est)>1)
        {
          DfD<-diag.level(Y,Phi.est,cov.shrinkage=cov.shrinkage)
        }else
        {
          DfD<-1
        }
      }
    }
    if(stop.crt[1]=="DfD")
    {
      nD<-1
      if(score.return)
      {
        score<-array(NA,c(n,max(nV),nD))
        dimnames(score)[[2]]<-paste0("visit",1:max(nV))
        dimnames(score)[[3]]<-paste0("D",1:nD)
        
        score[,,1]<-re1$score
        score.tmp<-score
      }
      
      DfD.tmp<-1
      while(DfD.tmp<DfD.thred)
      {
        re.tmp<-NULL
        try(tm.tmp<-system.time(re.tmp<-cap_Dk(Y,X,Phi0=Phi.est,method="CAP",OC=OC,cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,score.return=score.return,trace=FALSE,
                                               gamma0.mat=gamma0.mat,ninitial=ninitial,seed=seed)))
        if(is.null(re.tmp)==FALSE)
        {
          nD<-nD+1
          
          DfD<-diag.level(Y,cbind(Phi.est,re.tmp$gamma),cov.shrinkage=cov.shrinkage)
          DfD.tmp<-DfD$avg.level[nD]
          
          if(DfD.tmp<DfD.thred)
          {
            Phi.est<-cbind(Phi.est,re.tmp$gamma)
            beta.est<-cbind(beta.est,re.tmp$beta)
            beta0.random.est<-cbind(beta0.random.est,re.tmp$beta0.random)
            beta0.sigma2.est<-cbind(beta0.sigma2.est,re.tmp$beta0.sigma2)
            
            cp.time<-cbind(cp.time,as.numeric(tm.tmp[1:3]))
            
            if(verbose)
            {
              print(paste0("Component ",ncol(Phi.est)))
            }
            
            if(cov.shrinkage)
            {
              sk.out<-rbind(sk.out,re.tmp$shrinkage)
            }
            if(score.return)
            {
              score<-array(NA,c(n,max(nV),nD))
              dimnames(score)[[2]]<-paste0("visit",1:max(nV))
              dimnames(score)[[3]]<-paste0("D",1:nD)
              
              score[,,1:(nD-1)]<-score.tmp
              score[,,nD]<-re.tmp$score
              
              score.tmp<-score
            }
          }
        }else
        {
          break
        }
      }
      
      colnames(Phi.est)=colnames(beta.est)=colnames(beta0.random.est)=colnames(beta0.sigma2.est)<-paste0("D",1:ncol(Phi.est))
      rownames(Phi.est)<-paste0("V",1:p)
      rownames(beta.est)<-colnames(X[[1]])
      cp.time<-cbind(cp.time,apply(cp.time,1,sum))
      colnames(cp.time)<-c(paste0("D",1:ncol(Phi.est)),"Total")
      if(cov.shrinkage)
      {
        rownames(sk.out)<-paste0("D",1:ncol(Phi.est))
      }
      
      if(ncol(Phi.est)>1)
      {
        DfD<-diag.level(Y,Phi.est,cov.shrinkage=cov.shrinkage)
      }else
      {
        DfD<-1
      }
    }
    
    if(score.return)
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=Phi.est,beta=beta.est,beta0.random=beta0.random.est,beta0.sigma2=beta0.sigma2.est,
                 orthogonality=t(Phi.est)%*%Phi.est,DfD=DfD,shrinkage=sk.out,score=score,time=cp.time)
      }else
      {
        re<-list(gamma=Phi.est,beta=beta.est,beta0.random=beta0.random.est,beta0.sigma2=beta0.sigma2.est,
                 orthogonality=t(Phi.est)%*%Phi.est,DfD=DfD,score=score,time=cp.time)
      }
    }else
    {
      if(cov.shrinkage)
      {
        re<-list(gamma=Phi.est,beta=beta.est,beta0.random=beta0.random.est,beta0.sigma2=beta0.sigma2.est,
                 orthogonality=t(Phi.est)%*%Phi.est,DfD=DfD,shrinkage=sk.out,time=cp.time)
      }else
      {
        re<-list(gamma=Phi.est,beta=beta.est,beta0.random=beta0.random.est,beta0.sigma2=beta0.sigma2.est,
                 orthogonality=t(Phi.est)%*%Phi.est,DfD=DfD,time=cp.time)
      }
    }
    
    return(re)
  }
}
#################################################

#################################################
# Inference of beta
# based on bootstrap
cap.cov_beta_boot<-function(Y.cov,X,Y.n,gamma=NULL,boot=TRUE,sims=1000,boot.ci.type=c("perc","boot.se","bca"),conf.level=0.95,boot.seed=100,verbose=TRUE,max.itr=1000,tol=1e-4)
{
  # Y.cov: a list of covariance matrix of Y
  # X: covariates (longitudinal)
  # Y.n: n by nV matrix of sample size
  
  #--------------------------------------------
  n<-length(Y.cov)
  nV<-sapply(X,nrow)
  p<-dim(Y.cov[[1]])[1]
  
  q<-ncol(X[[1]])
  
  nT<-Y.n
  for(i in 1:n)
  {
    if(is.null(colnames(X[[i]]))==TRUE)
    {
      colnames(X[[i]])<-paste0("X",0:(q-1))
    }
  }
  #--------------------------------------------
  
  if(boot)
  {
    if(is.null(gamma))
    {
      stop("Error! Need gamma value.")
    }else
    {
      beta.boot<-matrix(NA,q,sims)
      beta0.random.boot<-matrix(NA,n,sims)
      beta0.sigma2.boot<-rep(NA,sims)
      
      for(b in 1:sims)
      {
        set.seed(boot.seed+b)
        # resample subjects with replacement
        idx.tmp<-sample(1:n,n,replace=TRUE)
        # print(idx.tmp)
        
        Ytmp.cov<-Y.cov[idx.tmp]
        Xtmp<-X[idx.tmp]
        nTtmp<-nT[idx.tmp,]
        
        re.tmp<-NULL
        try(re.tmp<-cap.cov_beta(Ytmp.cov,Xtmp,nTtmp,gamma,max.itr=max.itr,tol=tol,score.return=FALSE,trace=FALSE))
        if(is.null(re.tmp)==FALSE)
        {
          if(re.tmp$convergence==TRUE)
          {
            beta.boot[,b]<-re.tmp$beta 
            beta0.random.boot[,b]<-re.tmp$beta0.random
            beta0.sigma2.boot[b]<-re.tmp$beta0.sigma2
          }
        }
        
        if(verbose)
        {
          print(paste0("Bootstrap sample ",b))
        }
      }
      
      for(j in 1:q)
      {
        if(sum(is.na(beta.boot[j,]))>0)
        {
          itmp.nna<-which(is.na(beta.boot[j,])==FALSE)
          dis.cook<-cooks.distance(lm(beta.boot[j,]~1))
          cook.thred<-4/((length(itmp.nna)-2-2))
          itmp<-itmp.nna[which(dis.cook>cook.thred)]
          beta.boot[j,itmp]<-NA
          
          # print(length(itmp))
          
          while(length(itmp)>0)
          {
            itmp.nna<-which(is.na(beta.boot[j,])==FALSE)
            dis.cook<-cooks.distance(lm(beta.boot[j,]~1))
            cook.thred<-4/((length(itmp.nna)-2-2))
            itmp<-itmp.nna[which(dis.cook>cook.thred)]
            beta.boot[j,itmp]<-NA
            
            # print(length(itmp))
          }
        }
      }
      
      # summar
      beta.est<-apply(beta.boot,1,mean,na.rm=TRUE)
      beta.se<-apply(beta.boot,1,sd,na.rm=TRUE)
      beta.stat<-beta.est/beta.se
      beta.pv<-(1-pnorm(abs(beta.stat)))*2
      
      beta0.random.est<-apply(beta0.random.boot,1,mean,na.rm=TRUE)
      beta0.random.se<-apply(beta0.random.boot,1,sd,na.rm=TRUE)
      beta0.random.stat<-beta0.random.est/beta0.random.se
      beta0.random.pv<-(1-pnorm(abs(beta0.random.stat)))*2
      
      beta0.sigma2.est<-mean(beta0.sigma2.boot,na.rm=TRUE)
      beta0.sigma2.se<-sd(beta0.sigma2.boot,na.rm=TRUE)
      beta0.sigma2.stat<-beta0.sigma2.est/beta0.sigma2.se
      beta0.sigma2.pv<-(1-pnorm(abs(beta0.sigma2.stat)))*2
      
      if(sum(is.na(beta.boot))>0&boot.ci.type[1]=="bca")
      {
        boot.ci.type<-"perc"
      }
      if(boot.ci.type[1]=="bca")
      {
        beta.ci<-t(apply(beta.boot,1,BC.CI,sims=sims,conf.level=conf.level))
        beta0.random.ci<-t(apply(beta0.random.boot,1,BC.CI,sims=sims,conf.level=conf.level))
        beta0.sigma2.ci<-BC.CI(beta0.sigma2.boot,sims=sims,conf.level=conf.level)
      }
      if(boot.ci.type[1]=="perc")
      {
        beta.ci<-t(apply(beta.boot,1,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE))
        beta0.random.ci<-t(apply(beta0.random.boot,1,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE))
        beta0.sigma2.ci<-quantile(beta0.sigma2.boot,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      }
      if(boot.ci.type[1]=="boot.se")
      {
        zv<-qnorm(1-(1-conf.level)/2)
        beta.ci<-cbind(beta.est-zv*beta.se,beta.est+zv*beta.se)
        beta0.random.ci<-cbind(beta0.random.est-zv*beta0.random.se,beta0.random.est+zv*beta0.random.se)
        beta0.sigma2.ci<-c(beta0.sigma2.est-zv*beta0.sigma2.se,beta0.sigma2.est+zv*beta0.sigma2.se)
      }
      
      re.beta<-data.frame(Estimate=beta.est,SE=beta.se,statistic=beta.stat,pvalue=beta.pv,LB=beta.ci[,1],UB=beta.ci[,2])
      rownames(re.beta)<-colnames(X[[1]])
      
      re.beta0.random<-data.frame(Estimate=beta0.random.est,SE=beta0.random.se,statistic=beta0.random.stat,pvalue=beta0.random.pv,LB=beta0.random.ci[,1],UB=beta0.random.ci[,2])
      
      re.beta0.sigma2<-data.frame(Estimate=beta0.sigma2.est,SE=beta0.sigma2.se,statistic=beta0.sigma2.stat,pvalue=beta0.sigma2.pv,LB=beta0.sigma2.ci[1],UB=beta0.sigma2.ci[2])
      rownames(re.beta0.sigma2)<-"beta0.sigma2"
      
      re.inf<-list(beta=re.beta,beta0.random=re.beta0.random,beta0.sigma2=re.beta0.sigma2)
      re.boot<-list(beta=beta.boot,beta0.random=beta0.random.boot,beta0.sigma2=beta0.sigma2.boot)
      re<-list(Inference=re.inf,boot=re.boot)
      
      return(re)
    }
  }else
  {
    stop("Error!")
  }
}
cap_beta_boot<-function(Y,X,gamma=NULL,cov.shrinkage=TRUE,boot=TRUE,sims=1000,boot.ci.type=c("perc","boot.se","bca"),conf.level=0.95,boot.seed=100,verbose=TRUE,max.itr=1000,tol=1e-4)
{
  # Y: outcome list (longitudinal)
  # X: covariates (longitudinal)
  # gamma: linear projection
  
  #--------------------------------------------
  n<-length(Y)
  nV<-sapply(Y,length)
  p<-ncol(Y[[1]][[1]])
  
  q<-ncol(X[[1]])
  
  nT<-matrix(NA,n,max(nV))
  for(i in 1:n)
  {
    for(v in 1:nV[i])
    {
      nT[i,v]<-nrow(Y[[i]][[v]])
    }
    
    if(is.null(colnames(X[[i]]))==TRUE)
    {
      colnames(X[[i]])<-paste0("X",0:(q-1))
    }
  }
  
  if(min(nT,na.rm=TRUE)-5<p)
  {
    cov.shrinkage<-TRUE
  }
  #--------------------------------------------
  
  if(boot)
  {
    if(is.null(gamma))
    {
      stop("Error! Need gamma value.")
    }else
    {
      beta.boot<-matrix(NA,q,sims)
      beta0.random.boot<-matrix(NA,n,sims)
      beta0.sigma2.boot<-rep(NA,sims)
      
      for(b in 1:sims)
      {
        set.seed(boot.seed+b)
        # resample subjects with replacement
        idx.tmp<-sample(1:n,n,replace=TRUE)
        # print(idx.tmp)
        
        Ytmp<-Y[idx.tmp]
        Xtmp<-X[idx.tmp]
        
        re.tmp<-NULL
        try(re.tmp<-cap_beta(Ytmp,Xtmp,gamma,cov.shrinkage=cov.shrinkage,max.itr=max.itr,tol=tol,score.return=FALSE,trace=FALSE))
        if(is.null(re.tmp)==FALSE)
        {
          if(re.tmp$convergence==TRUE)
          {
            beta.boot[,b]<-re.tmp$beta 
            beta0.random.boot[,b]<-re.tmp$beta0.random
            beta0.sigma2.boot[b]<-re.tmp$beta0.sigma2
          }
        }
        
        if(verbose)
        {
          print(paste0("Bootstrap sample ",b))
        }
      }
      
      for(j in 1:q)
      {
        if(sum(is.na(beta.boot[j,]))>0)
        {
          itmp.nna<-which(is.na(beta.boot[j,])==FALSE)
          dis.cook<-cooks.distance(lm(beta.boot[j,]~1))
          cook.thred<-4/((length(itmp.nna)-2-2))
          itmp<-itmp.nna[which(dis.cook>cook.thred)]
          beta.boot[j,itmp]<-NA
          
          # print(length(itmp))
          
          while(length(itmp)>0)
          {
            itmp.nna<-which(is.na(beta.boot[j,])==FALSE)
            dis.cook<-cooks.distance(lm(beta.boot[j,]~1))
            cook.thred<-4/((length(itmp.nna)-2-2))
            itmp<-itmp.nna[which(dis.cook>cook.thred)]
            beta.boot[j,itmp]<-NA
            
            # print(length(itmp))
          }
        }
      }
      
      # summar
      beta.est<-apply(beta.boot,1,mean,na.rm=TRUE)
      beta.se<-apply(beta.boot,1,sd,na.rm=TRUE)
      beta.stat<-beta.est/beta.se
      beta.pv<-(1-pnorm(abs(beta.stat)))*2
      
      beta0.random.est<-apply(beta0.random.boot,1,mean,na.rm=TRUE)
      beta0.random.se<-apply(beta0.random.boot,1,sd,na.rm=TRUE)
      beta0.random.stat<-beta0.random.est/beta0.random.se
      beta0.random.pv<-(1-pnorm(abs(beta0.random.stat)))*2
      
      beta0.sigma2.est<-mean(beta0.sigma2.boot,na.rm=TRUE)
      beta0.sigma2.se<-sd(beta0.sigma2.boot,na.rm=TRUE)
      beta0.sigma2.stat<-beta0.sigma2.est/beta0.sigma2.se
      beta0.sigma2.pv<-(1-pnorm(abs(beta0.sigma2.stat)))*2
      
      if(sum(is.na(beta.boot))>0&boot.ci.type[1]=="bca")
      {
        boot.ci.type<-"perc"
      }
      if(boot.ci.type[1]=="bca")
      {
        beta.ci<-t(apply(beta.boot,1,BC.CI,sims=sims,conf.level=conf.level))
        beta0.random.ci<-t(apply(beta0.random.boot,1,BC.CI,sims=sims,conf.level=conf.level))
        beta0.sigma2.ci<-BC.CI(beta0.sigma2.boot,sims=sims,conf.level=conf.level)
      }
      if(boot.ci.type[1]=="perc")
      {
        beta.ci<-t(apply(beta.boot,1,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE))
        beta0.random.ci<-t(apply(beta0.random.boot,1,quantile,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE))
        beta0.sigma2.ci<-quantile(beta0.sigma2.boot,probs=c((1-conf.level)/2,1-(1-conf.level)/2),na.rm=TRUE)
      }
      if(boot.ci.type[1]=="boot.se")
      {
        zv<-qnorm(1-(1-conf.level)/2)
        beta.ci<-cbind(beta.est-zv*beta.se,beta.est+zv*beta.se)
        beta0.random.ci<-cbind(beta0.random.est-zv*beta0.random.se,beta0.random.est+zv*beta0.random.se)
        beta0.sigma2.ci<-c(beta0.sigma2.est-zv*beta0.sigma2.se,beta0.sigma2.est+zv*beta0.sigma2.se)
      }
      
      re.beta<-data.frame(Estimate=beta.est,SE=beta.se,statistic=beta.stat,pvalue=beta.pv,LB=beta.ci[,1],UB=beta.ci[,2])
      rownames(re.beta)<-colnames(X[[1]])
      
      re.beta0.random<-data.frame(Estimate=beta0.random.est,SE=beta0.random.se,statistic=beta0.random.stat,pvalue=beta0.random.pv,LB=beta0.random.ci[,1],UB=beta0.random.ci[,2])
      
      re.beta0.sigma2<-data.frame(Estimate=beta0.sigma2.est,SE=beta0.sigma2.se,statistic=beta0.sigma2.stat,pvalue=beta0.sigma2.pv,LB=beta0.sigma2.ci[1],UB=beta0.sigma2.ci[2])
      rownames(re.beta0.sigma2)<-"beta0.sigma2"
      
      re.inf<-list(beta=re.beta,beta0.random=re.beta0.random,beta0.sigma2=re.beta0.sigma2)
      re.boot<-list(beta=beta.boot,beta0.random=beta0.random.boot,beta0.sigma2=beta0.sigma2.boot)
      re<-list(Inference=re.inf,boot=re.boot)
      
      return(re)
    }
  }else
  {
    stop("Error!")
  }
}
#################################################



