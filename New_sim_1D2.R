setwd("~/Downloads")

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
library(ggpubr)
library(gridExtra)
# library(rstanarm)
# library(drc)
# library(sicegar)

can <- read.csv("HMS_LINCS_RPPA_Data_Normalized_(Innoscan_Mapix)_SRD_Sep_21.csv")
AZ_628=can[can$Small.Molecule.Name=='AZ-628',c(1:4,5:6,8,10,12,14:15,17,19,21:25)]
#AZ_628[,-(1:4)]<-scale(AZ_628[,-(1:4)])
d<-c(unique(AZ_628$Drug))

response <- read.csv("HMS_LINCS_Viability_Data_Normalized_SRD_Sep_21.csv")
AZ_628_res=response[response$Small.Molecule.Name=='AZ-628',]
#AZ_628_res[,-(1:4)]=scale(AZ_628_res[,-(1:4)])

#plot(1:7,mu0[5,],t="l")

X1<-apply(design(d=0.0032),2,function(x) (x-min(x))/(max(x)-min(x)))
X2<-apply(design(d=0.0100),2,function(x) (x-min(x))/(max(x)-min(x)))
X3<-apply(design(d=0.0316),2,function(x) (x-min(x))/(max(x)-min(x)))
X4<-apply(design(d=0.1000),2,function(x) (x-min(x))/(max(x)-min(x)))
X5<-apply(design(d=0.3160),2,function(x) (x-min(x))/(max(x)-min(x)))
X6<-apply(design(d=1.0000),2,function(x) (x-min(x))/(max(x)-min(x)))
X7<-apply(design(d=3.1600),2,function(x) (x-min(x))/(max(x)-min(x)))
#apply(X1,2,function(x) (x-min(x))/(max(x)-min(x)))

dta<-list(X1,X2,X3,X4,X5,X6,X7)
#dta=xfun(dta)
# true Values
a=.2
b=.8
c=2
theta=.6


mu01<-0.2+((0.8-0.2)/(1+(2/d)^0.6))
mu0<-matrix(0,nrow = 10,7)
set.seed(5)
for(i in 1:7) mu0[,i]<-msm::rtnorm(10,mu01[i],0.005,lower = 0)

#matplot(1:7,t(mu0),lty = 1,col=1:10,type = "l",lwd = 1.5,xlab="d")

beta<-seq(.15,.18,length.out = 14)#seq(.32,.65,length.out = 14)
beta_g<-seq(.16,.20,length.out = 14)#seq(.21,.62,length.out = 14)
# names(beta)<-c("" ,"b.1" , "b.2",  "b.3" , "b.4" , "b.5" , "b.6",  "b.7",  "b.8",  "b.9",  "b.10","b.11", "b.12", "b.13", "b.14")
# names(beta_g)<-c("" ,"b.1" , "b.2",  "b.3" , "b.4" , "b.5" , "b.6",  "b.7",  "b.8",  "b.9",  "b.10","b.11", "b.12", "b.13", "b.14")

al_int<-c()
al<-rep(0,10)
for (i in 1:7){
  al<-seq(1.51,1.51,length.out = 10)
  al_int<-cbind(al_int,al*i^1.2)
}
i=2

#exp(1/i^(20))

#matplot(1:7,t(al_int),lty = 1,col=1:50,type = "l",lwd = 1.5)
ga_int<-c()
ga<-rep(0,10)
for (i in 1:7){
  ga<-seq(1.31,1.61,length.out = 10)
  ga_int<-cbind(ga_int,ga*i^1.2)
}

#ga_int<-Ga


t<-1:3#c(10,24,48)
mu_new<-array(0,c(3,10,7))#t,i,d
for (k in 1:3){
  ui<-c()
  for (i in 1:7){
    ui<-cbind(ui,rnorm(10,mu0[,i]*exp(ga_int[,i]*(1-exp(-al_int[,i]*t[k]))),0.05))
  }
  mu_new[k,,]<-ui
}

AZ_628_res=response[response$Small.Molecule.Name=='AZ-628',]
t1=(unique(AZ_628_res$Time.Point..hr.))
y<-array(0,c(3,10,7))#t,i,d
for (j in 1:3){
  X00<-matrix(AZ_628_res[(AZ_628_res$Time.Point..hr. ==t1[j]),6],nrow = 10,byrow = T)#
  y[j,,]<-X00
}



#currently interpolating
y[1,7,]<-(y[1,8,]+y[1,6,])/2
y[1,10,]<-(y[1,9,]+y[1,1,])/2
y[3,9,]<-(y[3,8,]+y[3,10,])/2


mu_int<-mu_new
#al_new<-al_int
#ga_new<-ga_int

acc0<-0
acc1<-0
acc_alpha<-0
acc_gamma<-0
M=10000

sig1=0.05
sig2=0.05
sig0=0.05
sig1s<-c()
sig2s<-c()
sig0s<-c()

y=Y

# betaA1=betaA2=betaA3=betaA4=betaA5=betaA6=betaA7=betaS1<-matrix(0,14,M)
# betaG1=betaG2=betaG3=betaG4=betaG5=betaG6=betaG7=betaSG1<-matrix(0,14,M)

tempag1=tempag2=tempag3=tempag4=tempag5=tempag6=tempag7<-matrix(0.1,10,2)
temp_a1=matrix(al_int[,1],10,M)
temp_a2=matrix(al_int[,2],10,M)
temp_a3=matrix(al_int[,3],10,M)
temp_a4=matrix(al_int[,4],10,M)
temp_a5=matrix(al_int[,5],10,M)
temp_a6=matrix(al_int[,6],10,M)
temp_a7=matrix(al_int[,7],10,M)
temp_g1=matrix(ga_int[,1],10,M)
temp_g2=matrix(ga_int[,2],10,M)
temp_g3=matrix(ga_int[,3],10,M)
temp_g4=matrix(ga_int[,4],10,M)
temp_g5=matrix(ga_int[,5],10,M)
temp_g6=matrix(ga_int[,6],10,M)
temp_g7=matrix(ga_int[,7],10,M)
#mu_int=MUT

a1<-dta[[1]]%*%(beta)
a2<-dta[[2]]%*%(beta)
a3<-dta[[3]]%*%(beta)
a4<-dta[[4]]%*%(beta)
a5<-dta[[5]]%*%(beta)
a6<-dta[[6]]%*%(beta)
a7<-dta[[7]]%*%(beta)

g1<-dta[[1]]%*%(beta_g)
g2<-dta[[2]]%*%(beta_g)
g3<-dta[[3]]%*%(beta_g)
g4<-dta[[4]]%*%(beta_g)
g5<-dta[[5]]%*%(beta_g)
g6<-dta[[6]]%*%(beta_g)
g7<-dta[[7]]%*%(beta_g)


Xbind<-rbind(X1,X2,X3,X4,X5,X6,X7)
AAbind<-c()
GGbind<-c()
BetaA=c()
BetaG=c()

a=.1
b=.5
c=5
theta=2
#mu0

As<-c()
Bs<-c()
Cs<-c()
thetas<-c()
p=2

for (p in 2:M){
  print(p)
  
  a1 <- dta[[1]] %*% (beta)
  a2 <- dta[[2]] %*% (beta)
  a3 <- dta[[3]] %*% (beta)
  a4 <- dta[[4]] %*% (beta)
  a5 <- dta[[5]] %*% (beta)
  a6 <- dta[[6]] %*% (beta)
  a7 <- dta[[7]] %*% (beta)
  
  g1 <- dta[[1]] %*% (beta_g)
  g2 <- dta[[2]] %*% (beta_g)
  g3 <- dta[[3]] %*% (beta_g)
  g4 <- dta[[4]] %*% (beta_g)
  g5 <- dta[[5]] %*% (beta_g)
  g6 <- dta[[6]] %*% (beta_g)
  g7 <- dta[[7]] %*% (beta_g)
  
  
  BetaA = cbind(BetaA, beta)
  BetaG = cbind(BetaG, beta_g)
  
  threA <- rep(1, 10)
  threG <- rep(1, 10)
  threM <- rep(1, 10)
  pMu0<-priorMUpdate(a,b,c,theta,d[1:7])
  
  tempag1 <-alpha_gamma_update(y= y[,,1],mean_u0 = mu0[,1],priorA = a1,priorG = g1,alpha = temp_a1[,(p-1)],gamma = temp_g1[,(p-1)],thA=threA,thG=threG,sigma1 = sqrt(sig2),priorA1=a1,priorG1=g1)
  #mu_int[,,1] <-mutdi_update(y = y[,,1],mu = mu0[,1],al = tempag1[,1],ga = tempag1[,2],sigma1 = sqrt(sig1),sigma2 = sqrt(sig2))
  mu0[,1]<-mu0Update(y=y[,,1],mu0=mu0[,1],alpha=tempag1[,1],gamma=tempag1[,2],sigma2 = sqrt(sig2),pMu0=pMu0[1],pSigma0=sig0)
  #mu0[,1]<-mu0_update(mutdi = mu_int[,,1],al = tempag1[,1],ga = tempag1[,2],sigma1 = sqrt(sig1),sigma0 = sqrt(sig0),a=.5,b=1,c=3,theta=0.5,do=d[1])
  temp_a1[,p]=threA=tempag1[,1]
  temp_g1[,p]=threG=tempag1[,2]
  #threM<-mu0[,1]
  
  print(p)
  tempag2<-alpha_gamma_update(y=y[,,2],mean_u0 = mu0[,2],priorA = a2,priorG = g2,alpha = temp_a2[,(p-1)],gamma = temp_g2[,(p-1)],thA=threA,thG=threG,sigma1 = sqrt(sig2),priorA1=a1,priorG1=g1)
  mu0[,2]<-mu0Update(y=y[,,2],mu0=mu0[,2],alpha=tempag2[,1],gamma=tempag2[,2],sigma2 = sqrt(sig2),pMu0=pMu0[2],pSigma0=sig0)
  temp_a2[,p]=threA=tempag2[,1]
  temp_g2[,p]=threG=tempag2[,2]
  #threM<-mu0[,2]
  
  #print(p)
  tempag3<-alpha_gamma_update(y=y[,,3],mean_u0 = mu0[,3],priorA = a3,priorG = g3,alpha = temp_a3[,(p-1)],gamma = temp_g3[,(p-1)],thA=threA,thG=threG,sigma1 = sqrt(sig2),priorA1=a2,priorG1=g2)
  mu0[,3]<-mu0Update(y=y[,,3],mu0=mu0[,3],alpha=tempag3[,1],gamma=tempag3[,2],sigma2 = sqrt(sig2),pMu0=pMu0[3],pSigma0=sig0)
  temp_a3[,p]=threA=tempag3[,1]
  temp_g3[,p]=threG=tempag3[,2]
  # #threM<-mu0[,3]
  # 
  # #print(p)
  tempag4<-alpha_gamma_update(y=y[,,4],mean_u0 = mu0[,4],priorA = a4,priorG = g4,alpha = temp_a4[,(p-1)],gamma = temp_g4[,(p-1)],thA=threA,thG=threG,sigma1 = sqrt(sig2),priorA1=a3,priorG1=g3)
  mu0[,4]<-mu0Update(y=y[,,4],mu0=mu0[,4],alpha=tempag4[,1],gamma=tempag4[,2],sigma2 = sqrt(sig2),pMu0=pMu0[4],pSigma0=sig0)
  temp_a4[,p]=threA=tempag4[,1]
  temp_g4[,p]=threG=tempag4[,2]
  # #threM<-mu0[,4]
  # 
  # #print(p)
  tempag5<-alpha_gamma_update(y=y[,,5],mean_u0 = mu0[,5],priorA = a5,priorG = g5,alpha = temp_a5[,(p-1)],gamma = temp_g5[,(p-1)],thA=threA,thG=threG,sigma1 = sqrt(sig2),priorA1=a4,priorG1=g4)
  mu0[,5]<-mu0Update(y=y[,,5],mu0=mu0[,5],alpha=tempag5[,1],gamma=tempag5[,2],sigma2 = sqrt(sig2),pMu0=pMu0[5],pSigma0=sig0)
  temp_a5[,p]=threA=tempag5[,1]
  temp_g5[,p]=threG=tempag5[,2]
  # #threM<-mu0[,5]
  # 
  # #print(p)
  tempag6<-alpha_gamma_update(y=y[,,6],mean_u0 = mu0[,6],priorA = a6,priorG = g6,alpha = temp_a6[,(p-1)],gamma = temp_g6[,(p-1)],thA=threA,thG=threG,sigma1 = sqrt(sig2),priorA1=a5,priorG1=g5)
  mu0[,6]<-mu0Update(y=y[,,6],mu0=mu0[,6],alpha=tempag6[,1],gamma=tempag6[,2],sigma2 = sqrt(sig2),pMu0=pMu0[6],pSigma0=sig0)
  temp_a6[,p]=threA=tempag6[,1]
  temp_g6[,p]=threG=tempag6[,2]
  # #threM<-mu0[,6]
  # 
  # #print(p)
  tempag7<-alpha_gamma_update(y=y[,,7],mean_u0 = mu0[,7],priorA = a7,priorG = g7,alpha = temp_a7[,(p-1)],gamma = temp_g7[,(p-1)],thA=threA,thG=threG,sigma1 = sqrt(sig2),priorA1=a6,priorG1=g6)
  mu0[,7]<-mu0Update(y=y[,,7],mu0=mu0[,7],alpha=tempag7[,1],gamma=tempag7[,2],sigma2 = sqrt(sig2),pMu0=pMu0[7],pSigma0=sig0)
  temp_a7[,p]=threA=tempag7[,1]
  temp_g7[,p]=threG=tempag7[,2]
  # #threM<-mu0[,7]
  # 
  Al_pha<-cbind(temp_a1[,p],temp_a2[,p],temp_a3[,p],temp_a4[,p],temp_a5[,p],temp_a6[,p],temp_a7[,p])
  Ga_mma<-cbind(temp_g1[,p],temp_g2[,p],temp_g3[,p],temp_g4[,p],temp_g5[,p],temp_g6[,p],temp_g7[,p])
  
  #Mu0Sigma0<-updatepMu0Sigma0(mu0)
  Mu0Sigma0<-updateABCTheta(mu0,sig0,a,b,c,theta,d)
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
  sig2<-sigma0Update(y,mu0,Ga_mma,Al_pha,sig2)#sqrt(rinvgamma(1,shape =(1+210/2),scale =sig2_shape+1 ))
  sig2s<-c(sig2s,sig2)
  # 
  # 
  Gb<-matrix(cbind(log(temp_g1[,p]-1),
                   log(temp_g2[,p]-temp_g1[,p]),
                   log(temp_g3[,p]-temp_g2[,p]),
                   log(temp_g4[,p]-temp_g3[,p]),
                   log(temp_g5[,p]-temp_g4[,p]),
                   log(temp_g6[,p]-temp_g5[,p]),
                   log(temp_g7[,p]-temp_g6[,p])),nrow = 70,byrow = F)
  
  
  Ab<-matrix(cbind(log(temp_a1[,p]-1),
                   log(temp_a2[,p]-temp_a1[,p]),
                   log(temp_a3[,p]-temp_a2[,p]),
                   log(temp_a4[,p]-temp_a3[,p]),
                   log(temp_a5[,p]-temp_a4[,p]),
                   log(temp_a6[,p]-temp_a5[,p]),
                   log(temp_a7[,p]-temp_a6[,p])),nrow = 70,byrow = F)
  # 
  #  
  AAbind<-cbind(AAbind,Ab)
  GGbind<-cbind(GGbind,Gb)
  # 
  # #var(y-mean(y))
  # 
  capture.output(baa<-blasso(Xbind,Ab,RJ=F,lambda2 = 0,), file='NUL')#RJ=T
  beta<-colMeans(baa$beta[-(1:500),])#c(mean(baa$mu),colMeans(baa$beta[-(1:500),]))
  
  capture.output(bag<-blasso(Xbind,Gb,RJ=F,lambda2 = 0), file='NUL')#RJ=T
  beta_g<-colMeans(bag$beta[-(1:500),])#c(mean(bag$mu),colMeans(bag$beta[-(1:500),]))
}


trplot(temp_g1)
trplot(temp_g2)
trplot(temp_g3)
trplot(temp_g4)
trplot(temp_g5)
trplot(temp_g6)
trplot(temp_g7)

trplot(temp_a1)
trplot(temp_a2)
trplot(temp_a3)
trplot(temp_a4)
trplot(temp_a5)
trplot(temp_a6)
trplot(temp_a7)

AT

rowMeans(temp_a1[,-(1:M/2)])
rowMeans(temp_a2[,-(1:M/2)])
rowMeans(temp_a3[,-(1:M/2)])
rowMeans(temp_a4[,-(1:M/2)])
rowMeans(temp_a5[,-(1:M/2)])
rowMeans(temp_a6[,-(1:M/2)])
rowMeans(temp_a7[,-(1:M/2)])

rowMeans(temp_g1[,-(1:M/2)])
rowMeans(temp_g2[,-(1:M/2)])
rowMeans(temp_g3[,-(1:M/2)])
rowMeans(temp_g4[,-(1:M/2)])
rowMeans(temp_g5[,-(1:M/2)])
rowMeans(temp_g6[,-(1:M/2)])
rowMeans(temp_g7[,-(1:M/2)])

#geweke.plot(as.mcmc(temp_g1[,-(1:M/2)]))
geweke.diag(as.mcmc(t(temp_a1[,-(1:M/2)])))
geweke.diag(as.mcmc(t(temp_g1[,-(1:M/2)])))
geweke.diag(as.mcmc(As[-(1:M/2)]))
geweke.diag(as.mcmc(Bs[-(1:M/2)]))
geweke.diag(as.mcmc(Cs[-(1:M/2)]))
geweke.diag(as.mcmc(thetas[-(1:M/2)]))

par(mfrow=c(1,1))
plot(c(1:15),rowMeans(temp_a1[,-(1:M/2)]),typ = 'l',ylim=c(0,60),main = "Alpha1",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_a1[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_a1[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),AT[,1],pch=8)  

plot(c(1:15),rowMeans(temp_a2[,-(1:M/2)]),typ = 'l',ylim=c(0,60),main = "Alpha2",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_a2[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_a2[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),AT[,2],pch=8)

plot(c(1:15),rowMeans(temp_a3[,-(1:M/2)]),typ = 'l',ylim=c(0,60),main = "Alpha3",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_a3[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_a3[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),AT[,3],pch=8)

plot(c(1:15),rowMeans(temp_a4[,-(1:M/2)]),typ = 'l',ylim=c(0,60),main = "Alpha4",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_a4[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_a4[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),AT[,4],pch=8)

plot(c(1:15),rowMeans(temp_a5[,-(1:M/2)]),typ = 'l',ylim=c(0,60),main = "Alpha5",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_a5[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_a5[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),AT[,5],pch=8)

plot(c(1:15),rowMeans(temp_a6[,-(1:M/2)]),typ = 'l',ylim=c(0,60),main = "Alpha6",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_a6[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_a6[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),AT[,6],pch=8)

plot(c(1:15),rowMeans(temp_a7[,-(1:M/2)]),typ = 'l',ylim=c(0,60),main = "Alpha7",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_a7[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_a7[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),AT[,7],pch=8)


#gamma
par(mfrow=c(1,1))
plot(c(1:15),rowMeans(temp_g1[,-(1:M/2)]),typ = 'l',ylim=c(0,10),main = "Gamma1",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_g1[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_g1[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),GT[,1],pch=8)

plot(c(1:15),rowMeans(temp_g2[,-(1:M/2)]),typ = 'l',ylim=c(0,40),main = "Gamma2",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_g2[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_g2[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),GT[,2],pch=8)

plot(c(1:15),rowMeans(temp_g3[,-(1:M/2)]),typ = 'l',ylim=c(0,20),main = "Gamma3",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_g3[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_g3[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),GT[,3],pch=8)

plot(c(1:15),rowMeans(temp_g4[,-(1:M/2)]),typ = 'l',ylim=c(0,20),main = "Gamma4",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_g4[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_g4[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),GT[,4],pch=8)

plot(c(1:15),rowMeans(temp_g5[,-(1:M/2)]),typ = 'l',ylim=c(0,30),main = "Gamma5",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_g5[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_g5[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),GT[,5],pch=8)

plot(c(1:15),rowMeans(temp_g6[,-(1:M/2)]),typ = 'l',ylim=c(0,30),main = "Gamma6",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_g6[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_g6[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),GT[,6],pch=8)

plot(c(1:15),rowMeans(temp_g7[,-(1:M/2)]),typ = 'l',ylim=c(0,30),main = "Gamma7",ylab="Posterior Mean",xlab="Subjects")
points(c(1:15),apply(temp_g7[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:15),apply(temp_g7[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:15),GT[,7],pch=8)





quantile(temp_a1[1,-(1:M/5)],c(0.025,0.975))
# quantile(temp_a1[2,-(1:M/5)],c(0.025,0.975))
# quantile(temp_a1[3,-(1:M/5)],c(0.025,0.975))
# quantile(temp_a1[4,-(1:M/5)],c(0.025,0.975))
# quantile(temp_a1[5,-(1:M/5)],c(0.025,0.975))
# quantile(temp_a1[6,-(1:M/5)],c(0.025,0.975))
# quantile(temp_a1[7,-(1:M/5)],c(0.025,0.975))
# AT[,1]


quantile(sig0s[(2000:4999)],0.025)

# 
#dev.off()
par(mfrow=c(1,1))
traceplot(as.mcmc(sig0s[-(1:M/2)]),main="Sigma0")
#traceplot(as.mcmc(sig1s[-(1:1500)]),main="Sigma1")
traceplot(as.mcmc(sig2s[-(1:M/2)]),main="Sigma2")
traceplot(as.mcmc(As[-(1:M/2)]),main="A")
traceplot(as.mcmc(Bs[-(1:M/2)]),main="B")
traceplot(as.mcmc(Cs[-(1:M/2)]),main="C")
traceplot(as.mcmc(thetas[-(1:M/2)]),main="Theta")



mix_beta2ql<-function(x=non2) return( apply(x,2,function(x) (ifelse(length(x[x==0])>550,0,quantile(x[x!=0],0.025)))))         
mix_beta2qu<-function(x=non2) return( apply(x,2,function(x) (ifelse(length(x[x==0])>550,0,quantile(x[x!=0],0.975)))))         

log(rowMeans(temp_a1[,-(M/5)]))
log(AT[,1])


# MCMC
set.seed(111)
fit1<-cv.glmnet(x=X1,y=log(rowMeans(temp_a1[,-(M/5)])),
                alpha=.5,family="gaussian",nfolds = 3)
round(coef(fit1,s="lambda.min"),3)
beA[,1]

set.seed(118)
fit2<-cv.glmnet(x=X2,y=log(rowMeans(temp_a2[,-(M/5)])-rowMeans(temp_a1[,-(1:M/5)])),
                alpha=.5,family="gaussian",nfolds = 3)
round(coef(fit2,s="lambda.min"),3)
beA[,2]

set.seed(120)
fit3<-cv.glmnet(x=X3,y=log(rowMeans(temp_a3[,-(M/5)])-rowMeans(temp_a2[,-(1:M/5)])),
                alpha=.5,family="gaussian",nfolds = 3)
round(coef(fit3,s="lambda.min"),3)
beA[,3]

set.seed(126)
fit4<-cv.glmnet(x=X4,y=log(rowMeans(temp_a4[,-(M/5)])-rowMeans(temp_a3[,-(1:M/5)])),
                alpha=.5,family="gaussian",nfolds = 3)
round(coef(fit4,s="lambda.min"),3)
beA[,4]

set.seed(135)
fit5<-cv.glmnet(x=X5,y=log(rowMeans(temp_a5[,-(M/5)])-rowMeans(temp_a4[,-(1:M/5)])),
                alpha=.5,family="gaussian",nfolds = 3)
round(coef(fit5,s="lambda.min"),3)
beA[,5]

set.seed(154)
fit6<-cv.glmnet(X6,log(rowMeans(temp_a6[,-(1:M/5)])-rowMeans(temp_a5[,-(1:M/5)])), 
          alpha=0.5,family="gaussian",nfolds = 3)
round(coef(fit6,s="lambda.min"),3)
beA[,6]

set.seed(158)
fit7<-cv.glmnet(X7,log(rowMeans(temp_a7[,-(1:M/5)])-rowMeans(temp_a6[,-(1:M/5)])), 
                alpha=0.5,family="gaussian",nfolds = 3)
round(coef(fit7,s="lambda.min"),3)
beA[,7]


##################### gamma

set.seed(15)
fit11<-cv.glmnet(x=X1,y=log(rowMeans(temp_g1[,-(M/5)])),
                alpha=.5,family="gaussian",nfolds = 3)
round(coef(fit11,s="lambda.min"),3)
beG[,1]

#log(GT[,2]-GT[,1])
set.seed(32)
fit22<-cv.glmnet(x=X2,y=log(rowMeans(temp_g2[,-(M/5)])-rowMeans(temp_g1[,-(1:M/5)])),
                alpha=.5,family="gaussian",nfolds = 3)
round(coef(fit22,s="lambda.min"),3)
beG[,2]

set.seed(19)
fit33<-cv.glmnet(x=X3,y=log(rowMeans(temp_g3[,-(M/5)])-rowMeans(temp_g2[,-(1:M/5)])),
                alpha=.5,family="gaussian",nfolds = 3)
round(coef(fit33,s="lambda.min"),3)
beG[,3]

set.seed(23)
fit44<-cv.glmnet(x=X4,y=log(rowMeans(temp_g4[,-(M/5)])-rowMeans(temp_g3[,-(1:M/5)])),
                alpha=.5,family="gaussian",nfolds = 3)
round(coef(fit44,s="lambda.min"),3)
beG[,4]

set.seed(13)
fit55<-cv.glmnet(x=X5,y=log(rowMeans(temp_g5[,-(M/5)])-rowMeans(temp_g4[,-(1:M/5)])),
                alpha=.5,family="gaussian",nfolds = 3)
round(coef(fit55,s="lambda.min"),3)
beG[,5]

set.seed(11)
fit66<-cv.glmnet(X6,log(rowMeans(temp_g6[,-(1:M/5)])-rowMeans(temp_g5[,-(1:M/5)])), 
                alpha=0.5,family="gaussian",nfolds = 3)
round(coef(fit66,s="lambda.min"),3)
beG[,6]

set.seed(21)
fit77<-cv.glmnet(X7,log(rowMeans(temp_g7[,-(1:M/5)])-rowMeans(temp_g6[,-(1:M/5)])), 
                alpha=0.5,family="gaussian",nfolds = 3)
round(coef(fit77,s="lambda.min"),3)
beG[,7]

# A1bind<-matrix(cbind(log(AT[,1]),
#                     log(AT[,2]-AT[,1]),
#                     log(AT[,3]-AT[,2]),
#                     log(AT[,4]-AT[,3]),
#                     log(AT[,5]-AT[,4]),
#                     log(AT[,6]-AT[,5]),
#                     log(AT[,7]-AT[,6])),nrow = 70,byrow = F)
# Xbind<-rbind(X1,X2,X3,X4,X5,X6,X7)
# Gbind<-matrix(cbind(log(rowMeans(temp_g1[,-(1:9000)])),
#                     log(rowMeans(temp_g2[,-(1:9000)])-rowMeans(temp_g1[,-(1:9000)])),
#                     log(rowMeans(temp_g3[,-(1:9000)])-rowMeans(temp_g2[,-(1:9000)])),
#                     log(rowMeans(temp_g4[,-(1:9000)])-rowMeans(temp_g3[,-(1:9000)])),
#                     log(rowMeans(temp_g5[,-(1:9000)])-rowMeans(temp_g4[,-(1:9000)])),
#                     log(rowMeans(temp_g6[,-(1:9000)])-rowMeans(temp_g5[,-(1:9000)])),
#                     log(rowMeans(temp_g7[,-(1:9000)])-rowMeans(temp_g6[,-(1:9000)]))),nrow = 70,byrow = F)
# 
# 
# Abind<-matrix(cbind(log(rowMeans(temp_a1[,-(1:9000)])),
#                     log(rowMeans(temp_a2[,-(1:9000)])-rowMeans(temp_a1[,-(1:9000)])),
#                     log(rowMeans(temp_a3[,-(1:9000)])-rowMeans(temp_a2[,-(1:9000)])),
#                     log(rowMeans(temp_a4[,-(1:9000)])-rowMeans(temp_a3[,-(1:9000)])),
#                     log(rowMeans(temp_a5[,-(1:9000)])-rowMeans(temp_a4[,-(1:9000)])),
#                     log(rowMeans(temp_a6[,-(1:9000)])-rowMeans(temp_a5[,-(1:9000)])),
#                     log(rowMeans(temp_a7[,-(1:9000)])-rowMeans(temp_a6[,-(1:9000)]))),nrow = 70,byrow = F)




dim(Xbind)
dim(Abind)
dim(Gbind)

#### glmnet for beta(alpha)
set.seed(5)
fit<-cv.glmnet(Xbind,Abind,alpha=0.5,family="gaussian",nfolds = 3)
round(coef(fit,s="lambda.min"),3)

#### glmnet for beta(gamma)
set.seed(1)
fit0<-cv.glmnet(Xbind,Gbind,alpha=0.5,family="gaussian",nfolds = 3)
round(coef(fit0,s="lambda.min"),3)

mix_beta3ql<-function(x=non2,m,q) return( round (apply(x,2,function(x) (ifelse(length(x[x==0])>m,0,quantile(x[x!=0],q)))),3))         
mix_beta3qu<-function(x=non2,m,q) return( round(apply(x,2,function(x) (ifelse(length(x[x==0])>m,0,quantile(x[x!=0],q)))),3))         

#function to throw off non significant variables
mix_beta3<-function(x=non2,m) return( round(apply(x,2,function(x) (ifelse(length(x[x==0])>m,0,mean(x[x!=0])))),3))        


#################################################################################
#################### Posterior mean and Credible interval. ######################
#################################################################################

#ALPHA
set.seed(134)
capture.output(ba<-blasso(Xbind,Abind,T = 1000), file='NUL')

#betas 
non<-(ba$beta)

# Number of zero entries for each variables
table(unlist(apply(non,1,function(x) which(x==0))))

#Letting varible zero if there is more than 480 zeros

mix_beta3(non,m=471) #Posterior Mean
mix_beta3ql(non,m=471,q=0.025) # Lower credible interval
mix_beta3qu(non,m=471,q=0.975) # Upper credible interval

par(mfrow=c(1,1))
plot(c(1:14),mix_beta3(non,m=471),typ = 'l',ylim=c(-.1,.5),main = "Beta Alpha CI",ylab="Posterior Mean",xlab="Variables")
points(c(1:14),mix_beta3ql(non,m=471,q=0.025))
points(c(1:14),mix_beta3qu(non,m=471,q=0.975))
points(c(1:14),beA[-15,1],pch=8)


####   GAMMA
set.seed(49)
capture.output(ba0<-blasso(Xbind,Gbind,T = 5000), file='NUL')

#betas 
non0<-(ba0$beta)[-(1:4000),]

# Number of zero entries for each variables
table(unlist(apply(non0,1,function(x) which(x==0))))

#Letting varible zero if there is more than 500 zeros
mix_beta3(non0,m=490) #Posterior Mean
mix_beta3ql(non0,m=490,q=0.025) # Lower credible interval
mix_beta3qu(non0,m=490,q=0.975) # Upper credible interval

par(mfrow=c(1,1))

plot(c(1:14),mix_beta3(non0,m=490),typ = 'l',ylim=c(-1,1),main = "Beta Gamma CI",ylab="Posterior Mean",xlab="Variables")
points(c(1:14),mix_beta3ql(non0,m=490,q=0.025))
points(c(1:14),mix_beta3qu(non0,m=490,q=0.975))
points(c(1:14),rowMeans(beG[-15,]),pch=8)


#################################################################################
#############################  MAP and HPD interval #############################
#################################################################################


#ALPHA
set.seed(111)
capture.output(ba<-blasso(Xbind,Abind,T = 5000,RJ=F,verb=0), file='NUL')

non<-(ba$beta)[-(1:2000),]
hpda<-HPDinterval(as.mcmc(non),0.95)# HPD interval
estimate_mode(non)# MAP estimate
hpda[,1]# Lower limit of HPD interval
hpda[,2]# Upper limit of HPD interval

par(mfrow=c(1,1))
plot(c(1:14),colMode(non),typ = 'l',ylim=c(-.1,.5),main = "Beta Alpha HPD",ylab="Posterior Mean",xlab="Variables")
points(c(1:14),hpda[,1])
points(c(1:14),hpda[,2])
points(c(1:14),rowMeans(beA[-15,]),pch=8)

####   GAMMA
set.seed(222)
#betas GAMMA
capture.output(ba0<-blasso(Xbind,Gbind,T = 5000,RJ=F), file='NUL')
non0<-(ba0$beta)[-(1:1000),]

hpdg<-HPDinterval(as.mcmc(non0),0.95)# HPD interval
estimate_mode(non0)# MAP estimate
hpdg[,1]# Lower limit of HPD interval
hpdg[,2]# Upper limit of HPD interval

par(mfrow=c(1,1))
plot(c(1:14),colMode(non0),typ = 'l',ylim=c(-1,1),main = "Beta Gamma HPD",ylab="Posterior Mean",xlab="Variables")
points(c(1:14),hpdg[,1])
points(c(1:14),hpdg[,2])
points(c(1:14),rowMeans(beG[-15,]),pch=8)

#### Postrior Mean of beta(alpha)
rowMeans(BetaA[,-(1:M/2)])

#### HPD interval
hpd<-HPDinterval(as.mcmc(t(BetaA[,-(1:M/2)])))
par(mfrow=c(1,1))
plot(c(1:14),rowMeans(BetaA[,-(1:M/2)]),typ = 'l',ylim=c(-5,5),main = "Beta Alpha CI",ylab="Posterior Mean",xlab="Variables")
points(c(1:14),apply(BetaA[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:14),apply(BetaA[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:14),rowMeans(beA),pch=8)

### Postrior Mean of beta(gamma)
rowMeans(BetaG[,-(1:M/2)])

#### HPD interval
hpd1<-HPDinterval(as.mcmc(t(BetaG[,-(1:M/2)])))
par(mfrow=c(1,1))
plot(c(1:14),rowMeans(BetaG[,-(1:M/2)]),typ = 'l',ylim=c(-5,5),main = "Beta Gamma CI",ylab="Posterior Mean",xlab="Variables")
points(c(1:14),apply(BetaG[,-(1:M/2)],1,function (x) quantile(x,0.025)))
points(c(1:14),apply(BetaG[,-(1:M/2)],1,function (x) quantile(x,0.975)))
points(c(1:14),rowMeans(beG),pch=8)

# 
# trplot2(BetaA)
# trplot2(BetaG)
# 
# 
# 
# # Correlation panel
# panel.cor <- function(x, y){
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   r <- round(cor(x, y), digits=2)
#   txt <- paste0("R = ", r)
#   cex.cor <- 0.8/strwidth(txt)
#   text(0.5, 0.5, txt,col = ifelse(as.double(substr(txt,5,8))>0.70,'red','black'),cex = 1.5)
# }
# # Customize upper panel
# upper.panel<-function(x, y){
#   points(x,y, pch = 19)
# }
# # Create the plots
# pairs(Xbind, 
#       lower.panel = panel.cor,
#       upper.panel = upper.panel)
# tex
# ifelse(as.double(substr(tex,5,8))>0.40,'red','black')
# 
# 
# 
# z1=seq(1,0,length.out = 10000)
# y1=z1^z1
# 
# plot(z1,y1)
# 
# install.packages("rstanarm")
# install.packages("bayess")
# library(bayess)
# library(rstanarm)
# fit<-BayesReg(Ab,Xbind)
# fit$postmeancoeff
library(ggplot2)
pl=data.frame(x=1:15,y=rowMeans(temp_a1[,-(1:M/2)]),z=AT[,1],
              l=apply(temp_a1[,-(1:M/2)],1,function (x) quantile(x,0.025)),
              u=apply(temp_a1[,-(1:M/2)],1,function (x) quantile(x,0.975)))

gplot<-function(x,pre,low,upp,tru,name="Alpha 1"){
  pl=data.frame(x=x,y=pre,z=tru,
                l=low,
                u=upp)
  
  ggplot(pl) + 
    geom_line(aes(x,y,colour="blue")) +
    geom_point(aes(x,z,colour="darkred")) +
    theme_gray() +
    labs(title = paste0("Prediction including 95%CI ",name),x="Dose",y="Values") +
    geom_ribbon(aes(x,ymin=l, ymax=u,color="lightblue"),fill="lightblue", alpha=0.3)+
    scale_x_continuous(breaks=x,limit=c(1,15))+
    scale_color_manual(name=NULL,values = c("blue", "red","lightblue"),labels=c("Predicted", "True","95% CI"))+
    theme(plot.title = element_text(hjust=0.5),legend.position = c(0.95, 0.95),legend.justification = c("right", "top"),legend.background = element_rect(fill=NA))
  
  ggsave(paste0(name,".png"))
}
gplot(1:15,rowMeans(temp_a1[,-(1:M/2)]),apply(temp_a1[,-(1:M/2)],1,function (x) quantile(x,0.025)),
      apply(temp_a1[,-(1:M/2)],1,function (x) quantile(x,0.975)),AT[,1])

df=temp_a1[,-(1:M/2)]
tracegg<-function(df,variable="Alpha1"){
  n=dim(df)[1]
  m=dim(df)[2]
  plots <- vector("list", n)
  for (i in seq_along(plots))
    plots[[i]] <- ggplot(data.frame(df[i,]),aes(x=1:m,y=df[i,])) +
    geom_line()+
    labs(subtitle = paste0("Subject ",i,sep=""),x="Iteration",y="Values")+
    scale_x_continuous(breaks=c(1,m/2,m))+
    #scale_y_continuous(breaks=as.character(round(min(df[i,])-1,2)),as.character(round(max(df[i,])+1,2)),as.character(round(mean(df[i,]),2)))+
    theme(plot.subtitle = element_text(hjust=0.5,size=8),axis.title.x = element_text(size=8),
          axis.title.y = element_text(size=8),axis.text.x = element_text(size=5),
          axis.text.y = element_text(size=5))
  figure=do.call(ggarrange, plots)
  annotate_figure(figure,top = text_grob(paste0("Trace plot ",variable),face = "bold"))
  ggsave(paste0("Trace plot ",variable,".png"))
}
tracegg(df)
# figure=do.call(ggarrange, c(ncol = 2, nrow = 3,plots))
# 
# figure <- ggarrange(kp,kp,
#                     ncol = 1, nrow = 2)

# annotate_figure(figure,
#                 top = text_grob("Visualizing mpg", color = "red", face = "bold", size = 14),
#                 bottom = text_grob("Data source: \n mtcars data set", color = "blue",
#                                    hjust = 1, x = 1, face = "italic", size = 10),
#                 left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
#                 right = "I'm done, thanks :-)!",
#                 fig.lab = "Figure 1", fig.lab.face = "bold")



#scale_color_manual(name="Legend",values = c("blue", "red","green"),labels=c("True", "Predicted","95% CI")) +
  #t
# ggplot(by_year_percentage, aes(x=arrivaldate)) +
#   geom_line(aes(y=deathpercentage, color = "deathpercentage"), size = 1.5) +
#   geom_line(aes(y=tamponadepercentage, color = "tamponadepercentage"), size = 1.5) +
#   geom_line(aes(y=protaminepercentage, color = "protaminepercentage"), size = 1.5) +
#   scale_color_manual(name = "Group",
#                      values = c( "deathpercentage" = "blue", "tamponadepercentage" = "red", "protaminepercentage" = "orange"),
#                      labels = c("deathpercentage", "tamponadepercentage", "protaminepercentage"))
