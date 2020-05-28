#setwd("~/Downloads")

library(fda)
library(msm)
library(MASS)
library(mvtnorm)
library(coda)
library(truncnorm)
library(invgamma)
library(graphics)
library(EnvStats)
library(monomvn)
library(robustbase)
library(glmnet)
library(bayestestR)
library(truncdist)
library(bayess)
library(TruncatedNormal)
library(imputeTS)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(parallel)
library(MASS)
library(foreach)
library(doParallel)
# library(rstanarm)
# library(drc)
# library(sicegar)

can <- read.csv("HMS_LINCS_RPPA_Data_Normalized_(Innoscan_Mapix)_SRD_Sep_21.csv")
AZ_628=can[can$Small.Molecule.Name=='AZ-628',c(1:4,5:6,8,10,12,14:15,17,19,21:25)]
#AZ_628[,-(1:4)]<-scale(AZ_628[,-(1:4)])
d<-c(unique(AZ_628$Drug))
cell<-c(unique(AZ_628$Cell.Line.Name))

response <- read.csv("HMS_LINCS_Viability_Data_Normalized_SRD_Sep_21.csv")
AZ_628_res=response[response$Small.Molecule.Name=='AZ-628',]

# Initial Values
a=.3
b=.7
c=1
theta=.8

mu0p1<-0.3+((0.7-0.3)/(1+(1/d)^.8))
mu0p<-matrix(0,nrow = 10,7)
set.seed(5)
for(i in 1:7) mu0p[,i]<-msm::rtnorm(10,mu0p1[i],0.005,lower = 0,u=1)

matplot(1:7,t(mu0p),lty = 1,col=1:10,type = "l",lwd = 1.5,xlab="d")
mu0p

x1<-designCell(data=AZ_628,cell=cell[1])
x2<-designCell(data=AZ_628,cell=cell[2])
x3<-designCell(data=AZ_628,cell=cell[3])
x4<-designCell(data=AZ_628,cell=cell[4])
x5<-designCell(data=AZ_628,cell=cell[5])
x6<-designCell(data=AZ_628,cell=cell[6])
x7<-designCell(data=AZ_628,cell=cell[7])
x8<-designCell(data=AZ_628,cell=cell[8])
x9<-designCell(data=AZ_628,cell=cell[9])
x10<-designCell(data=AZ_628,cell=cell[10])

#dta<-list(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10)

beta<-seq(.15,.18,length.out = 14)
beta_g<-seq(.16,.20,length.out = 14)

al_int<-c()
al<-rep(0,10)
for (i in 1:7){
  al<-seq(1.51,1.61,length.out = 10)
  al_int<-cbind(al_int,al*i^1.2)
}

ga_int<-c()
ga<-rep(0,10)
for (i in 1:7){
  ga<-seq(1.31,1.61,length.out = 10)
  ga_int<-cbind(ga_int,ga*i^1.2)
}

M=10000

sig2=0.05
sig0=0.05
sig2s<-c()
sig0s<-c()

y=Y

resAlpha<-array(0,c(M,7,10))
resGamma<-array(0,c(M,7,10))
resMu0<-array(0,c(M,7,10))
resAlpha[1,,]<-t(al_int)
resGamma[1,,]<-t(ga_int)
resMu0[1,,]<-t(mu0p)

AAbind<-c()
GGbind<-c()
BetaA<-c()
BetaG<-c()

As<-c()
Bs<-c()
Cs<-c()
thetas<-c()
mu0S<-array(0,c(M,10,7))
mu0S[1,,]<-mu0

s2g=0.5
s2a=0.5

numCores <- detectCores()
numCores
Xbind1<-rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
Xbind<-apply(Xbind1,2,function(x) (x-min(x))/(max(x)-min(x)))

dta=lapply(seq_along(seq(1,70,by =7)), function(i) Xbind[(i:(i+6)),])

i=2
#p=2
#p=p+1
for (p in 2:M){
  print(p)
  
  XAbeta<-matrix(Xbind%*%beta,nrow = 10,byrow = T)
  XGbeta<-matrix(Xbind%*%beta_g,nrow = 10,byrow = T)
  
  BetaA = cbind(BetaA, beta)
  BetaG = cbind(BetaG, beta_g)
  
  pMu0<-muMu0(a,b,c,theta,d)
  
  registerDoParallel(numCores)
  Results=foreach(i=1:10) %dopar% {
    posAlpGamMu0(y=y[2:3,i,],y0=y[1,i,],mu0=resMu0[p-1,,i],alpha=resAlpha[p-1,,i],gamma=resGamma[p-1,,i],XAbeta=XAbeta[i,],
               XGbeta=XGbeta[i,],sigma=sqrt(sig2),betaASD=s2a,betaGSD=s2g,pMu0,pSigma0=sqrt(sig0))
  }
  
  
  
  resAlpha[p,,]<-do.call(cbind, lapply(Results, function(x) x[, 1]))
  resGamma[p,,]<-do.call(cbind, lapply(Results, function(x) x[, 2]))
  resMu0[p,,]<-do.call(cbind, lapply(Results, function(x) x[, 3]))
  
  #Al_pha<-cbind(temp_a1[,p],temp_a2[,p],temp_a3[,p],temp_a4[,p],temp_a5[,p],temp_a6[,p],temp_a7[,p])
  #Ga_mma<-cbind(temp_g1[,p],temp_g2[,p],temp_g3[,p],temp_g4[,p],temp_g5[,p],temp_g6[,p],temp_g7[,p])
  
  sig2<-posSigma2(y,t(resMu0[p,,]),t(resAlpha[p,,]),t(resGamma[p,,]),sig2)
  sig2s<-c(sig2s,sig2)
  
  Mu0Sigma0<-posABCTheSig0(t(resMu0[p,,]),sig0,a,b,c,theta)
  theta<-Mu0Sigma0[4]
  a<-Mu0Sigma0[1]
  b<-Mu0Sigma0[2]
  c<-Mu0Sigma0[3]
  sig0<-Mu0Sigma0[5]
  
  As<-c(As,a)
  Bs<-c(Bs,b)
  Cs<-c(Cs,c)
  thetas<-c(thetas,theta)
  sig0s<-c(sig0s,sig0)
  
  #sig2_shape<-sum((mu_int-y)^2)/2
  #print(p)
  
  #                                      
  # 
  Ab=foreach(i=2:7,.combine = c) %dopar%{
    log(Results2[[1]][i,]-Results2[[1]][i-1,])
  }
  Ab<-c(log(Results2[[1]][1,]-1),Ab)
  
  Gb=foreach(i=2:7,.combine = c) %dopar%{
    log(Results2[[2]][i,]-Results2[[2]][i-1,])
  }
  Gb<-c(log(Results2[[2]][1,]-1),Gb)
  
  # 
  #  
  AAbind<-cbind(AAbind,Ab)
  GGbind<-cbind(GGbind,Gb)
  # 
  # #var(y-mean(y))
  # 
  capture.output(baa<-blasso(Xbind,Ab,RJ=F,lambda2 = 0), file='NUL')#RJ=T
  beta<-colMeans(baa$beta[-(1:500),])#c(mean(baa$mu),colMeans(baa$beta[-(1:500),]))
  s2a=(sqrt(mean(baa$s2)))
  
  capture.output(bag<-blasso(Xbind,Gb,RJ=F,lambda2 = 0), file='NUL')#RJ=T
  beta_g<-colMeans(bag$beta[-(1:500),])#c(mean(bag$mu),colMeans(bag$beta[-(1:500),]))
  s2g=(sqrt(mean(bag$s2)))
  
  #stopImplicitCluster()
}

# mu0m<-matrix(0,10,7)
# for (s in M/2:(M-1)){
#   mu0m<-mu0m+mu0S[s,,]
# }
# mu0m<-mu0m/(M-1-(M/2))

# Amean=cbind(rowMeans(temp_a1[,-(1:M/2)]),rowMeans(temp_a2[,-(1:M/2)]),rowMeans(temp_a3[,-(1:M/2)]),
#       rowMeans(temp_a4[,-(1:M/2)]),rowMeans(temp_a5[,-(1:M/2)]),rowMeans(temp_a6[,-(1:M/2)]),
#       rowMeans(temp_a7[,-(1:M/2)]))
# 
# Gmean=cbind(rowMeans(temp_g1[,-(1:M/2)]),rowMeans(temp_g2[,-(1:M/2)]),rowMeans(temp_g3[,-(1:M/2)]),
#       rowMeans(temp_g4[,-(1:M/2)]),rowMeans(temp_g5[,-(1:M/2)]),rowMeans(temp_g6[,-(1:M/2)]),
#       rowMeans(temp_g7[,-(1:M/2)]))
betaap<-rowMeans(BetaG[,-(1:M/2)])
betagp<-rowMeans(BetaA[,-(1:M/2)])
Sig<-mean(sqrt(sig2s[-(1:M/2)]))

ATHat=matrix(exp(dta5[[1]]%*%betaap)+1,nrow = 5,7)
for(i in 2:7){
  ATHat[,i]=ATHat[,(i-1)]+exp(dta5[[i]]%*%betaap)
}

GTHat=matrix(exp(dta5[[1]]%*%betagp)+1,nrow = 5,7)
for(i in 2:7){
  GTHat[,i]=GTHat[,(i-1)]+exp(dta5[[i]]%*%betagp)
}

y0p<-Z[1,11:15,]
r48<-Z[2,11:15,]
r72<-Z[3,11:15,]

nm48<-as.vector(y0p*GTHat^(1-ATHat^(-1)))
p48dist<-rtmvnorm(100000,mu =nm48,sigma = diag(Sig,35), lb = rep(0,35),ub = rep(1,35))
p48<-matrix(colMeans(p48dist),5,7,byrow = F)
p48LCI<-matrix(apply(p48dist,2,function (x) quantile(x,0.025)),5,7,byrow = F)
p48UCI<-matrix(apply(p48dist,2,function (x) quantile(x,0.975)),5,7,byrow = F)

nm72<-as.vector(y0p*GTHat^(1-ATHat^(-2)))
p72dist<-rtmvnorm(100000,mu =nm72,sigma = diag(Sig,35), lb = rep(0,35),ub = rep(1,35))
p72<-matrix(colMeans(p72dist),5,7,byrow = F)
p72LCI<-matrix(apply(p72dist,2,function (x) quantile(x,0.025)),5,7,byrow = F)
p72UCI<-matrix(apply(p72dist,2,function (x) quantile(x,0.975)),5,7,byrow = F)

MAE48<-mean(abs(r48-p48))
MAE72<-mean(abs(r72-p72))

TMAE<-sum(MAE48,MAE72)/2


PMSE48<-(mean((r48-p48)^2))
PMSE72<-(mean((r72-p72)^2))

PRMSE<-(sum(PMSE48,PMSE72)/2)
                
#currently interpolating
# y[1,7,]<-(y[1,8,]+y[1,6,])/2
# y[1,10,]<-(y[1,9,]+y[1,1,])/2
# y[3,9,]<-(y[3,8,]+y[3,10,])/2

