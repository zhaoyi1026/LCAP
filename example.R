###########################################
# longitudinal covariance matrix regression

# R example
###########################################

library("mvtnorm")

rm(list=ls())

###################################
# parameter setting

p<-20
q<-3

# l2-norm 1
set.seed(100)
gamma.mat0<-matrix(runif(p),nrow=p,ncol=p)

gamma.mat<-qr.Q(qr(gamma.mat0))
for(j in 1:p)
{
  if(gamma.mat[which.max(abs(gamma.mat[,j])),j]<0)
  {
    gamma.mat[,j]<-(-gamma.mat[,j])
  }
}
Gamma<-gamma.mat
# t(Gamma)%*%Gamma
# Gamma%*%t(Gamma)

beta1.base<-seq(3,-1,length.out=5)
beta1.vec<-c(beta1.base,seq(-1.5,-3,length.out=p-length(beta1.base)))

beta2.base<-c(0,-0.5,0,0.5,0)
beta2.vec<-c(beta2.base,rep(0,p-length(beta2.base)))
beta3.base<-c(0,0.5,0,-0.25,0)
beta3.vec<-c(beta3.base,rep(0,p-length(beta3.base)))

beta.mat<-rbind(beta1.vec,beta2.vec,beta3.vec)
rownames(beta.mat)<-paste0("beta",1:q)
colnames(beta.mat)<-paste0("dim",1:p)

# beta2=0 corresponding component beta SD
beta.nr.sd<-0.1

# beta2 \neq 0 corresponding component SD of beta in longitudinal
beta.rv.sd<-0.1
###################################

###################################
# generate data

n<-100
nV.m<-5
nT.m<-100

# # of visits of each subject
set.seed(100)
nV<-round(rnorm(n,nV.m,sd=1))
# # of observations of each subject at each visit
set.seed(100)
nT<-matrix(NA,n,max(nV))     
for(i in 1:n)
{
  nT[i,1:nV[i]]<-round(rnorm(nV[i],mean=nT.m,sd=5))
}

# generate beta0
beta0.mat<-matrix(NA,n,ncol(beta.mat))
colnames(beta0.mat)<-colnames(beta.mat)
for(j in 1:ncol(beta.mat))
{
  if(beta.mat[2,j]==0)
  {
    set.seed(100)
    beta0.mat[,j]<-rnorm(n,mean=beta.mat[1,j],sd=beta.nr.sd)
  }else
  {
    set.seed(100)
    beta0.mat[,j]<-rnorm(n,mean=beta.mat[1,j],sd=beta.rv.sd)
  }
}

# generate data
set.seed(100)
X<-vector("list",length=n)
Y<-vector("list",length=n)
Sigma<-vector("list",length=n)
delta<-vector("list",length=n)
for(i in 1:n)
{
  # X1: binary and time-invariant
  # X2: continuous and time-varying
  X[[i]]<-cbind(rep(1,nV[i]),rep(rbinom(1,size=1,prob=0.5),nV[i]),rnorm(nV[i],mean=0,sd=0.5))
  rownames(X[[i]])<-paste0("visit",1:nV[i])
  colnames(X[[i]])<-paste0("X",1:q)
  
  Y[[i]]<-vector("list",length=nV[i])
  Sigma[[i]]<-array(NA,c(p,p,nV[i]))
  delta[[i]]<-matrix(NA,nV[i],p)
  for(v in 1:nV[i])
  {
    for(j in 1:ncol(beta.mat))
    {
      if(beta.mat[2,j]!=0)
      {
        delta[[i]][v,j]<-exp(t(X[[i]][v,])%*%c(beta0.mat[i,j],beta.mat[-1,j]))
      }else
      {
        delta[[i]][v,j]<-exp(t(X[[i]][v,])%*%c(beta0.mat[i,j],beta.mat[-1,j]))
      }
    }
    Sigma[[i]][,,v]<-Gamma%*%diag(delta[[i]][v,])%*%t(Gamma)
    
    Y[[i]][[v]]<-rmvnorm(n=nT[i,v],mean=rep(0,p),sigma=Sigma[[i]][,,v])
  }
}
###################################

###################################
# run method

source("Longitudinal_HDCAP.R")

# method parameters
cov.shrinkage<-TRUE
max.itr<-1000
tol<-1e-4
trace<-TRUE
score.return<-TRUE

stop.crt<-"DfD"
DfD.thred<-2
# nD<-5
verbose<-TRUE

boot<-TRUE
sims<-500
boot.ci.type<-c("boot.se")
conf.level<-0.95

re<-NULL
try(re<-capReg(Y,X,stop.crt="DfD",DfD.thred=DfD.thred,cov.shrinkage=TRUE,max.itr=max.itr,tol=tol,score.return=score.return,trace=trace,verbose=verbose))

if(is.null(re)==FALSE)
{
  re.boot<-vector("list",length=ncol(re$gamma))
  for(jj in 1:ncol(re$gamma))
  {
    re.boot[[jj]]<-NULL
    try(re.boot[[jj]]<-cap_beta_boot(Y,X,gamma=re$gamma[,jj],cov.shrinkage=TRUE,boot=boot,sims=sims,boot.ci.type=boot.ci.type,conf.level=conf.level,verbose=FALSE))
    
    print(paste0("Component ",jj))
  }
}
###################################
save.image("example.RData")

