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
library(nleqslv)
#if (dir.exists(paste0("/work/statgrads/krag57/DSim"))==F) dir.create(paste0("/work/statgrads/krag57/DSim"))
can <- read.csv("HMS_LINCS_RPPA_Data_Normalized_(Innoscan_Mapix)_SRD_Sep_21.csv")
AZ_628=can[can$Small.Molecule.Name=='AZ-628',c(1:4,5:6,8,10,12,14:15,17,19,21:25)]

source("Parallelized simulation.R")
# can <- read.csv("/home/statgrads/krag57/HMS_LINCS_RPPA_Data_Normalized_(Innoscan_Mapix)_SRD_Sep_21.csv")
# AZ_628=can[can$Small.Molecule.Name=='AZ-628',c(1:4,5:6,8,10,12,14:15,17,19,21:25)]

d<-c(unique(AZ_628$Drug))
cell<-c(unique(AZ_628$Cell.Line.Name))

response <- read.csv("HMS_LINCS_Viability_Data_Normalized_SRD_Sep_21.csv")
AZ_628_res=response[response$Small.Molecule.Name=='AZ-628' & response $Time.Point..hr==48,c(1,6)]
matrix(AZ_628_res[,2],ncol=7,byrow = T)


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

mu01<-0.1+((0.6-0.1)/(1+(2/d)^0.6))
mu015<-matrix(0,nrow = 15,7)
set.seed(5)
for(i in 1:7) mu015[,i]<-msm::rtnorm(15,mu01[i],0.005,lower = 0,u=1)

matplot(1:7,t(mu015),lty = 1,col=1:15,type = "l",lwd = 1.5,xlab="d")
mu0=mu015[1:10,]

Z=Y<-array(0,c(3,15,7))
# y0t<-c()
# for (i in 1:15){
#   y0i<-(msm::rtnorm(7,mu015[i,],0.005,lower = 0,u=1))
#   y0t=rbind(y0t,y0i)
# }
y0<-matrix(msm::rtnorm(105,mu015,0.005,lower = 0,u=1),nrow = 15,byrow = F)
Z[1,,]<-y0
for (k in 1:2){
  yt<-c()
  for (i in 1:15){
    yi<-(msm::rtnorm(7,y0[i,]*t(GT[,i])^(1-t(AT[,i])^(-k)),0.005,lower = 0,u=1))
    yt=rbind(yt,yi)
  }
  Z[k+1,,]<-yt
}
Z[2,,]  

Y=Z[,1:10,]

matplot(t(matrix(AZ_628_res[,2],ncol=7,byrow = T)),t="l")

# 1/as.vector((apply(Y,c(2,3),var)))
# plot(as.vector((apply(Y,c(2,3),var))))
# plot(density(1/as.vector((apply(Y,c(2,3),var)))))
# fitdistr(1/as.vector((apply(Y,c(2,3),var))),densfun = "gamma")
sd((apply(Y,c(2,3),var)))^2
mean((apply(Y,c(2,3),var)))/2
#a = 3.92901 and b = 0.150112
fn <- function(x) {
  rate <- x[2]/ (x[1]-1) - mean((apply(Y,c(2,3),var)))/2
  shape <- x[2]^2/((x[1]-1)^2*(x[1]-2))-sd((apply(Y,c(2,3),var))/2)^2
  return(c(rate, shape))
}
nleqslv(c(2.1,0.1), fn)$x
fitmodel <- nls(y ~a+(b-a)/(1 + (c/x)^o ),data = cbind.data.frame(y=Y[1,1,],x=d),start =list(a=0.1,b=0.4,c=2,o=1))
# summary(fitmodel)
# Formula: y ~ a + (b - a)/(1 + (c/x)^o)
# 
# Parameters:
#   Estimate Std. Error t value Pr(>|t|)   
# a  0.10320    0.01353   7.625  0.00468 **
# b  0.68431    0.26452   2.587  0.08129 . 
# c  3.04655    4.48353   0.679  0.54555   
# o  0.58286    0.16531   3.526  0.03875 * 
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.008085 on 3 degrees of freedom
# 
# Number of iterations to convergence: 6 
# Achieved convergence tolerance: 3.889e-06
mean(resid(fitmodel))
sd(resid(fitmodel))*sqrt(2)
#let shape is 3 and scale is =2*0.008085^2

#####################################  FANOVA   ######################################
######################################################################################
dosenames <- c("d1  ", "d2", "d3 ", "d4", "d5","d6","d7")
FA=Z[3,1:10,]
#  Set up a design matrix having a column for (the grand mean, and
#    a column for (dose. Add a dummy contraint
#    observation

dmat<-matrix(0,7,8)
dmat[,1] <- 1
dmat[1:7,2:8]<-diag(1,7,7)
#  labels for weather zones

dlabels <- vector("list",8)
dlabels[[1]] <- "Constant"
dlabels[[2]] <- "d1"
dlabels[[3]] <- "d2"
dlabels[[4]] <- "d3"
dlabels[[5]] <- "d4"
dlabels[[6]] <- "d5"
dlabels[[7]] <- "d6"
dlabels[[8]] <- "d7"

d36<-matrix(1,1,8)
d36[1]<-0
dmat   <- rbind(dmat, d36)

subjectrange <- c(1,10)
subjectbasis8 <- create.bspline.basis(subjectrange, 4)
smoothdList <- smooth.basis(1:10,FA,subjectbasis8, fdnames=list("Subject", "dose", "response"))
subjecttempfd  <- smoothdList$fd

coef   <- subjecttempfd$coefs
coef8 <- cbind(coef,matrix(0,4,1))
subjecttempfd$coefs <- coef8

pd <- 8
xfddlist <- vector("list",pd)
for (j in 1:pd) xfddlist[[j]] <- dmat[,j]
#  set up the basis for (the regression functions

betacbasis  <- create.constant.basis(subjectrange)
betadbasis  <- create.bspline.basis(subjectrange,4)

betadlist <- vector("list",pd)
betadlist[[1]] <- betacbasis
for (j in 2:pd) betadlist[[j]] <- betadbasis
#  compute regression coefficient functions and
#  predicted functions=
fRegressdList <- fRegress(subjecttempfd, xfddlist, betadlist)
#  plot regression functions
# 
# par(mfrow=c(3,3))
# for (j in 1:pd) {
#   betaestParfdj <- betadestlist[[j]]
#   plot(betaestParfdj$fd, xlab="Subject", ylab="response")
#   title(dlabels[[j]])
# }
#  set up predicted functions
yhatdfdobj <- fRegressdList$yhatfdobj

yhatdmat    <- eval.fd(1:10, yhatdfdobj$fd)
ymatd       <- eval.fd(1:10, subjecttempfd)
SSE <- sum((ymatd[,1:7] - yhatdmat[,1:7])^2)
SSE



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

s2g=1
s2a=1
S2A<-c(s2a)
S2G<-c(s2g)
numCores <- detectCores()
numCores

i=2
#p=2
#p=p+1

for (p in 2:M){
  print(p)
  #source("/home/statgrads/krag57/Parallelized simulation5.R")
  XAbeta<-matrix(Xbind%*%beta,nrow = 7,byrow = F)
  XGbeta<-matrix(Xbind%*%beta_g,nrow = 7,byrow = F)
  
  BetaA = cbind(BetaA, beta)
  BetaG = cbind(BetaG, beta_g)
  
  pMu0<-muMu0(a,b,c,theta,d)
  
  registerDoParallel(numCores)
  Results=foreach(i=1:10) %dopar% {
    pos_DAlpGamMu0(y=y[2:3,i,],y0=y[1,i,],mu0=resMu0[p-1,,i],dalpha=resDAlpha[p-1,,i],dgamma=resDGamma[p-1,,i],
                   XAbeta=XAbeta[,i],XGbeta=XGbeta[,i],sigma=sqrt(sig2),betaASD=s2a,betaGSD=s2g,pMu0,pSigma0=sqrt(sig0))
  }
  stopImplicitCluster()
  
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
  
  sig2<-posSigma2(y,mu0 = t(resMu0[p,,]),Alpha = t(resAlpha[p,,]),Gamma = t(resGamma[p,,]),Sigma2 = sig2)
  sig2s<-c(sig2s,sig2)
 
  Mu0Sigma0<-posABCTheSig0(mu0 = t(resMu0[p,,]),sig0 = sig0,a = a,b = b,c = c,theta = theta)
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
  
  capture.output(baa<-bridge(Xbind,matrix(log(resDAlpha[p,,]), ncol = 1),RJ=F,ab=c(3.93,0.150)), file='NUL')
  beta<-colMeans(baa$beta[-(1:500),])
  s2a=(sqrt(mean(baa$s2)))
  
  capture.output(bag<-bridge(Xbind,matrix(log(resDGamma[p,,]), ncol = 1),RJ=F,ab=c(3.93,0.150)), file='NUL')
  beta_g<-colMeans(bag$beta[-(1:500),])
  s2g=(sqrt(mean(bag$s2)))
  S2A<-c(S2A,s2a)
  S2G<-c(S2G,s2g)
}


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
tracePlotSub(resAlpha[(M/2):M,,],1)
# tracePlotSub(resAlpha[(M/2):M,,],2)
# tracePlotSub(resAlpha[(M/2):M,,],3)
# tracePlotSub(resAlpha[(M/2):M,,],4)
# tracePlotSub(resAlpha[(M/2):M,,],5)
# tracePlotSub(resAlpha[(M/2):M,,],6)
# tracePlotSub(resAlpha[(M/2):M,,],7)
# 
# 
tracePlotSub(resGamma[(M/2):M,,],1)
# tracePlotSub(resGamma[(M/2):M,,],2)
# tracePlotSub(resGamma[(M/2):M,,],3)
# tracePlotSub(resGamma[(M/2):M,,],4)
# tracePlotSub(resGamma[(M/2):M,,],5)
# tracePlotSub(resGamma[(M/2):M,,],6)
# tracePlotSub(resGamma[(M/2):M,,],7)
# 
# betaap<-rowMeans(BetaA[,-(1:M/2)])
# betagp<-rowMeans(BetaG[,-(1:M/2)])
# Sig<-(mean(sqrt(sig2s[-(1:M/2)])))
# 
par(mfrow=c(1,1))
traceplot(as.mcmc(sig0s[-(1:M/2)]),main="Sigma0")
quantile(sqrt(sig0s[-(1:M/2)]),c(0.025,0.975))
traceplot(as.mcmc(sig2s[-(1:M/2)]),main="Sigma2")
quantile((sig2s[-(1:M/2)]),c(0.025,0.975))
traceplot(as.mcmc(As[-(1:M/2)]),main="A")
quantile((As[-(1:M/2)]),c(0.025,0.975))
traceplot(as.mcmc(Bs[-(1:M/2)]),main="B")
quantile((Bs[-(1:M/2)]),c(0.025,0.975))
traceplot(as.mcmc(Cs[-(1:M/2)]),main="C")
quantile((Cs[-(1:M/2)]),c(0.025,0.975))
traceplot(as.mcmc(thetas[-(1:M/2)]),main="Theta")
quantile((thetas[-(1:M/2)]),c(0.025,0.975))
traceplot(as.mcmc(S2A[-(1:M/2)]),main="var1")
quantile((S2A[-(1:M/2)]),c(0.025,0.975))
traceplot(as.mcmc(S2G[-(1:M/2)]),main="var2")
quantile((S2G[-(1:M/2)]),c(0.025,0.975))
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


###########################################################################################
################################### Prediction ############################################
###########################################################################################

al=apply(resAlpha,c(2,3),function (x) quantile(x,0.025))
am=apply(resAlpha,c(2,3),mean)
at=AT[,1:10]
au=apply(resAlpha,c(2,3),function (x) quantile(x,0.975))


gplotpred3(t(am),t(al),t(au),t(at),variable="Subject")


gl=apply(resGamma,c(2,3),function (x) quantile(x,0.025))
gm=apply(resGamma,c(2,3),mean)
gt=GT[,1:10]
gu=apply(resGamma,c(2,3),function (x) quantile(x,0.975))
gplotpred3(t(gm),t(gl),t(gu),t(gt),variable="Gamma")

AlphaHat<-array(0,c(7,5,dim(BetaA[,-(1:M/2)])[2]))
for (r in 1:dim(BetaA[,-(1:M/2)])[2]){
  DA=matrix(0,nrow = 7,5)
  for(i in 1:5){
    DA[,i]=exp(dta15[[10+i]]%*%BetaA[,(M/2)+r])
  }
  
  AThat=matrix(0,nrow = 7,5)
  AThat[1,]<-DA[1,]+1
  for (i in 2:7){
    AThat[i,]<-AThat[i-1,]+DA[i,]
  }
  AlphaHat[,,r]<-AThat
}

GammaHat<-array(0,c(7,5,dim(BetaG[,-(1:M/2)])[2]))
for (r in 1:dim(BetaG[,-(1:M/2)])[2]){
  DG=matrix(0,nrow = 7,5)
  for(i in 1:5){
    DG[,i]=exp(dta15[[10+i]]%*%BetaG[,(M/2)+r])
  }
  
  GThat=matrix(0,nrow = 7,5)
  GThat[1,]<-DG[1,]+1
  for (j in 2:7){
    GThat[j,]<-GThat[j-1,]+DG[j,]
  }
  GammaHat[,,r]<-GThat
}


an=apply(AlphaHat, c(1,2), mean)
anl=apply(AlphaHat, c(1,2), function(x) quantile(x,0.025))
anu=apply(AlphaHat, c(1,2), function(x) quantile(x,0.975))
AT[,11:15]
gplotpred3(t(an),t(anl),t(anu),t(AT[,11:15]),variable="Alpha")


gn=apply(GammaHat, c(1,2), mean)
gnl=apply(GammaHat, c(1,2), function(x) quantile(x,0.025))
gnu=apply(GammaHat, c(1,2), function(x) quantile(x,0.975))
GT[,11:15]
gplotpred3(t(gn),t(gnl),t(gnu),t(GT[,11:15]),variable="Gamma")


ac=cbind(am,an)
acl=cbind(al,anl)
acu=cbind(au,anu)
AT
gplotpred3(t(ac),t(acl),t(acu),t(AT),variable="Alpha")

gc=cbind(gm,gn)
gcl=cbind(gl,gnl)
gcu=cbind(gu,gnu)
GT
gplotpred3(t(gc),t(gcl),t(gcu),t(GT),variable="Gamma")


y0p<-Z[1,11:15,]
r48<-Z[2,11:15,]
r72<-Z[3,11:15,]



y0pp<-matrix(as.vector(y0p),nrow = 7,ncol = 5,byrow = T)
pre48<-array(0,c(5,7,dim(GammaHat)[3]))
for (r in 1:dim(GammaHat)[3]){
  Sig<-sqrt(sig2s[(M/2)+r])
  # yt<-c()
  # for (i in 1:5){
  #   yi<-(msm::rtnorm(7,y0pp[,i]*t(GammaHat[,i,r])^(1-t(AlphaHat[,i,r])^(-1)),Sig,lower = 0,u=1))
  #   yt=rbind(yt,yi)
  # }
  # Z[k+1,,]<-yt
  p481dist<-rtmvnorm(1,mu =as.vector(y0pp*GammaHat[,,r]^(1-AlphaHat[,,r]^(-1))),sigma = diag(Sig,35), lb = rep(0,35),ub = rep(1,35))
  p481<-matrix(p481dist,5,7,byrow = T)
  pre48[,,r]<-p481
}

p48=apply(pre48, c(1,2), mean)
p48l=apply(pre48, c(1,2), function(x) quantile(x,0.025))
p48u=apply(pre48, c(1,2), function(x) quantile(x,0.975))
r48
gplotpred3(p48,p48l,p48u,r48,variable="Y @ t=48hr")

pre72<-array(0,c(5,7,dim(GammaHat)[3]))
for (r in 1:dim(GammaHat)[3]){
  Sig<-sqrt(sig2s[(M/2)+r])
  p721dist<-rtmvnorm(1,mu =as.vector(y0pp*GammaHat[,,r]^(1-AlphaHat[,,r]^(-2))),sigma = diag(Sig,35), lb = rep(0,35),ub = rep(1,35))
  p721<-matrix(p721dist,5,7,byrow = T)
  pre72[,,r]<-p721
}

p72=apply(pre72, c(1,2), mean)
p72l=apply(pre72, c(1,2), function(x) quantile(x,0.025))
p72u=apply(pre72, c(1,2), function(x) quantile(x,0.975))
r72
gplotpred3(p72,p72l,p72u,r72,variable="Y @ t=72hr")


sigmoid<-function(theta,x) 1/(1+exp(-theta*x))#-theta*x

cost<-function(x,theta,y) 0.5*(sigmoid(theta,x)-y)^2

x=rnorm(100, 0,1)
y=sample(c(0,1), 100,r=T)
theta=1
cost(x,theta,y)

J<-function(theta){
  x=rnorm(100, 0,1)
  y=sample(c(0,1), 100,r=T)
  return (mean(cost(x,theta,y)))
}
X=seq(10,50,1)
Y=sapply(X, J,simplify = TRUE, USE.NAMES = TRUE)
plot(X,Y,t="l")
