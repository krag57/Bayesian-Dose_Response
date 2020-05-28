setwd("~/Documents/GitHub/Bayesian-Dose_Response")

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
#if (dir.exists(paste0("/work/statgrads/krag57/DSim"))==F) dir.create(paste0("/work/statgrads/krag57/DSim"))
can <- read.csv("HMS_LINCS_RPPA_Data_Normalized_(Innoscan_Mapix)_SRD_Sep_21.csv")
AZ_628=can[can$Small.Molecule.Name=='AZ-628',c(1:4,5:6,8,10,12,14:15,17,19,21:25)]

source("Parallelized simulation.R")
# can <- read.csv("/home/statgrads/krag57/HMS_LINCS_RPPA_Data_Normalized_(Innoscan_Mapix)_SRD_Sep_21.csv")
# AZ_628=can[can$Small.Molecule.Name=='AZ-628',c(1:4,5:6,8,10,12,14:15,17,19,21:25)]

d<-c(unique(AZ_628$Drug))
cell<-c(unique(AZ_628$Cell.Line.Name))

# response <- read.csv("HMS_LINCS_Viability_Data_Normalized_SRD_Sep_21.csv")
# AZ_628_res=response[response$Small.Molecule.Name=='AZ-628',]

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

Xbind1<-rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
Xbind15<-apply(rbind(Xbind1,matrix(rtruncnorm(490,a=0,mean = 5,sd = 2),nrow = 35,ncol = 14)),
               2,function(x) (x-min(x))/(max(x)-min(x)))
Xbind<-Xbind15[1:70,]
dta15<-lapply(seq_along(seq(1,105,by =7)), function(i) Xbind15[(i:(i+6)),])
dta<-lapply(seq_along(seq(1,70,by =7)), function(i) Xbind[(i:(i+6)),])

###################################################################################################
###################################################################################################

betaa<-rep(0,14)
betaa[1:4]<-c(0.15,.25,.1,.2)

betag<-rep(0,14)
betag[1:4]<-c(0.15,.1,.2,.25)

DAT=matrix(0,nrow = 7,15)
for(i in 1:15){
  DAT[,i]=exp(dta15[[i]]%*%betaa)
}

DGT=matrix(0,nrow = 7,15)
for(i in 1:15){
  DGT[,i]=exp(dta15[[i]]%*%betag)
}

AT=matrix(0,nrow = 7,15)
AT[1,]<-DAT[1,]+1
for (i in 2:7){
  AT[i,]<-AT[i-1,]+DAT[i,]
}

GT=matrix(0,nrow = 7,15)
GT[1,]<-DGT[1,]+1
for (i in 2:7){
  GT[i,]<-GT[i-1,]+DGT[i,]
}

mu01<-0.2+((0.6-0.2)/(1+(2/d)^0.6))
mu015<-matrix(0,nrow = 15,7)
set.seed(5)
for(i in 1:7) mu015[,i]<-msm::rtnorm(15,mu01[i],0.005,lower = 0,u=1)

matplot(1:7,t(mu015),lty = 1,col=1:15,type = "l",lwd = 1.5,xlab="d")
mu0=mu015[1:10,]

Z=Y<-array(0,c(3,15,7))
y0<-matrix(msm::rtnorm(105,mu015,0.005,lower = 0,u=1),nrow = 15,byrow = F)
Z[1,,]<-y0
for (k in 1:2){
  yi<-matrix(msm::rtnorm(105,y0*t(GT)^(1-t(AT)^(-k)),0.005,lower = 0,u=1),nrow = 15,byrow = F)
  Z[k+1,,]<-yi
}
Y=Z[,1:10,]

########################################################################################################
########################################################################################################

beta<-seq(.15,.18,length.out = 14)
beta_g<-seq(.16,.20,length.out = 14)

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

beta<-seq(.15,.18,length.out = 14)
beta_g<-seq(.16,.20,length.out = 14)

al_int<-c()
for (i in 1:10){
  al_int<-cbind(al_int,sample(seq(1.10,1.51,length.out = 20),size = 7))
}

ga_int<-c()
for (i in 1:10){
  ga_int<-cbind(ga_int,sample(seq(1.10,1.51,length.out = 20),size = 7))
}

M=10000

sig2=0.05
sig0=0.05
sig2s<-c()
sig0s<-c()

y=Y

resDAlpha<-array(0,c(M,7,10))
resDGamma<-array(0,c(M,7,10))
resAlpha<-array(0,c(M,7,10))
resGamma<-array(0,c(M,7,10))
resMu0<-array(0,c(M,7,10))
resDAlpha[1,,]<-al_int
resDGamma[1,,]<-ga_int
resMu0[1,,]<-mu0p

AAbind<-c()
GGbind<-c()
BetaA<-c()    
BetaG<-c()

As<-c()
Bs<-c()
Cs<-c()
thetas<-c()

s2g=0.5
s2a=0.5
S2A<-c(s2a)
S2G<-c(s2g)
numCores <- detectCores()
numCores

i=2
#p=2
#p=p+1
registerDoParallel(numCores)
for (p in 2:M){
  print(p)
  #source("/home/statgrads/krag57/Parallelized simulation5.R")
  XAbeta<-matrix(Xbind%*%beta,nrow = 7,byrow = F)
  XGbeta<-matrix(Xbind%*%beta_g,nrow = 7,byrow = F)
  
  BetaA = cbind(BetaA, beta)
  BetaG = cbind(BetaG, beta_g)
  
  pMu0<-muMu0(a,b,c,theta,d)
  
  
  Results=foreach(i=1:10) %dopar% {
    pos_DAlpGamMu0(y=y[2:3,i,],y0=y[1,i,],mu0=resMu0[p-1,,i],dalpha=resDAlpha[p-1,,i],dgamma=resDGamma[p-1,,i],
                   XAbeta=XAbeta[,i],XGbeta=XGbeta[,i],sigma=sqrt(sig2),betaASD=s2a,betaGSD=s2g,pMu0,pSigma0=sqrt(sig0))
  }
    
  resDAlpha[p,,]<-do.call(cbind, lapply(Results, function(x) x[, 1]))
  resDGamma[p,,]<-do.call(cbind, lapply(Results, function(x) x[, 2]))
  resMu0[p,,]<-do.call(cbind, lapply(Results, function(x) x[, 3]))
  
  resAlpha[p,1,]<-resDAlpha[p,1,]+1
  for (i in 2:7){
    resAlpha[p,i,]<-resAlpha[p,i-1,]+resDAlpha[p,i,]
  }
  
  resGamma[p,1,]<-resDGamma[p,1,]+1
  for (i in 2:7){
    resGamma[p,i,]<-resGamma[p,i-1,]+resDGamma[p,i,]
  }
  
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
  
  AAbind<-cbind(AAbind,matrix(log(resDAlpha[p,,]), ncol = 1))
  GGbind<-cbind(GGbind,matrix(log(resDGamma[p,,]), ncol = 1))
  
  capture.output(baa<-bridge(Xbind,matrix(log(resDAlpha[p,,]), ncol = 1),RJ=F,ab=c(1, .05)), file='NUL')
  beta<-colMeans(baa$beta[-(1:500),])
  s2a=(sqrt(mean(baa$s2)))
  
  capture.output(bag<-bridge(Xbind,matrix(log(resDGamma[p,,]), ncol = 1),RJ=F,ab=c(1,.05)), file='NUL')
  beta_g<-colMeans(bag$beta[-(1:500),])
  s2g=(sqrt(mean(bag$s2)))
  S2A<-c(S2A,s2a)
  S2G<-c(S2G,s2g)
}
stopImplicitCluster()

apply(resAlpha,c(2,3),function (x) quantile(x,0.025))
apply(resAlpha,c(2,3),mean)
AT[,1:10]
apply(resAlpha,c(2,3),function (x) quantile(x,0.975))

apply(resGamma,c(2,3),function (x) quantile(x,0.025))
apply(resGamma,c(2,3),mean)
GT[,1:10]
apply(resGamma,c(2,3),function (x) quantile(x,0.975))

# apply(resDAlpha,c(2,3),mean)
# AT
# alphaGamma(resAlpha[(M/2):M,,],1,1)
# traceplot(as.mcmc(alphaGamma(resAlpha[(M/2):M,,],7,1)))
# betaap<-rowMeans(BetaA[,-(1:M/2)])
# betagp<-rowMeans(BetaG[,-(1:M/2)])
# Sig<-(mean(sqrt(sig2s[-(1:M/2)])))
# 
par(mfrow=c(1,1))
traceplot(as.mcmc(sig0s[-(1:M/2)]),main="Sigma0")
#traceplot(as.mcmc(sig1s[-(1:1500)]),main="Sigma1")
traceplot(as.mcmc(sig2s[-(1:M/2)]),main="Sigma2")
traceplot(as.mcmc(As[-(1:M/2)]),main="A")
traceplot(as.mcmc(Bs[-(1:M/2)]),main="B")
traceplot(as.mcmc(Cs[-(1:M/2)]),main="C")
traceplot(as.mcmc(thetas[-(1:M/2)]),main="Theta")
traceplot(as.mcmc(S2A[-(1:M/2)]),main="var1")
traceplot(as.mcmc(S2G[-(1:M/2)]),main="var2")

#### Postrior Mean of beta(alpha)
rowMeans(BetaA[,-(1:M/2)])

#### HPD interval
hpd<-HPDinterval(as.mcmc(t(BetaA[,-(1:M/2)])))
par(mfrow=c(1,1))
plot(c(1:14),rowMeans(BetaA[,-(1:M/2)]),typ = 'l',ylim=c(-5,5),main = "Beta Alpha CI",ylab="Posterior Mean",xlab="Variables")
points(c(1:14),apply(BetaA[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:14),apply(BetaA[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:14),(betaa),pch=8)

### Postrior Mean of beta(gamma)
rowMeans(BetaG[,-(1:M/2)])

#### HPD interval
hpd1<-HPDinterval(as.mcmc(t(BetaG[,-(1:M/2)])))
par(mfrow=c(1,1))
plot(c(1:14),rowMeans(BetaG[,-(1:M/2)]),typ = 'l',ylim=c(-5,5),main = "Beta Gamma CI",ylab="Posterior Mean",xlab="Variables")
points(c(1:14),apply(BetaG[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:14),apply(BetaG[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:14),betag,pch=8)

