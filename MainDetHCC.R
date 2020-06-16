#setwd("~/Documents/GitHub/Bayesian-Dose_Response")

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

M=10000
source("/home/statgrads/krag57/Function Det HCC.R")


can <- read.csv("/home/statgrads/krag57/HMS_LINCS_RPPA_Data_Normalized_(Innoscan_Mapix)_SRD_Sep_21.csv")
AZ_628=can[can$Small.Molecule.Name=='AZ-628',c(1:4,5:6,8,10,12,14:15,17,19,21:25)]
d<-c(unique(AZ_628$Drug))
cell<-c("C32","COLO 858","RVH-421","WM-115","WM1552C","LOXIMVI","MMAC-SF","MZ7-mel","K2","SKMEL28")

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
set.seed(77777)
Xbind15<-apply(rbind(Xbind1,matrix(rtruncnorm(490,a=0,mean = 5,sd = 2),nrow = 35,ncol = 14)),
               2,function(x) (x-min(x))/(max(x)-min(x)))
Xbind<-Xbind15[1:70,]

#registerDoParallel(10)
EntireData <- list()

for(h in 1:50){
  ###########################################################################################################
  ###########################################    Simulation   ###############################################
  ###########################################################################################################
  
  betaa<-rep(0,14)
  set.seed(h+1000)
  betaa[c(sample(1:14,size = 4+h%%4,replace = F))]<-sample(c(0.20,0.25,0.15,0.10),size = 4+h%%4,replace = T)
  
  betag<-rep(0,14)
  set.seed(h+2000)
  betag[c(sample(1:14,size = 4+h%%4,replace = F))]<-sample(c(0.20,0.25,0.15,0.10),size = 4+h%%4,replace = T)
  
  set.seed(h+1500)
  at<-runif(1,0,0.5)
  bt<-runif(1,at,.9)
  ct<-runif(1,0,3)
  thetat<-runif(1,0,2.5)
  
  DAT<-matrix(exp(Xbind15%*%betaa),nrow = 7)
  DGT<-matrix(exp(Xbind15%*%betag),nrow = 7)
  
  
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
  
  mu01<-at+((bt-at)/(1+(ct/d)^thetat))
  mu015<-matrix(rep(mu01,15),nrow = 15,7,byrow = T)
  
  mu0=mu015[1:10,]
  
  Z=Y<-array(0,c(3,15,7))
  # y0t<-c()
  # for (i in 1:15){
  #   y0i<-(msm::rtnorm(7,mu015[i,],0.005,lower = 0,u=1))
  #   y0t=rbind(y0t,y0i)
  # }
  set.seed(h+4000)
  y0<-matrix(msm::rtnorm(105,mu015,0.005,lower = 0,u=1),nrow = 15,byrow = F)
  Z[1,,]<-y0
  set.seed(h+4500)
  for (k in 1:2){
    yt<-c()
    for (i in 1:15){
      yi<-(msm::rtnorm(7,y0[i,]*t(GT[,i])^(1-t(AT[,i])^(-k)),0.005,lower = 0,u=1))
      yt=rbind(yt,yi)
    }
    Z[k+1,,]<-yt
  }
  
  Y=Z[,1:10,]
  
  cenAG<-mean((apply(Y,c(2,3),var)))/2
  
  ABCThetaTrue<-c(at,bt,ct,thetat,0.005,0.005,.6,0.6)
  
  ##########################  FANOVA   #######################
  ############################################################
  FA=Z[1,1:10,]
  dmat<-matrix(0,7,8)
  dmat[,1] <- 1
  dmat[1:7,2:8]<-diag(1,7,7)
  
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
  for (l in 1:pd) xfddlist[[l]] <- dmat[,l]
  #  set up the basis for (the regression functions
  
  betacbasis  <- create.constant.basis(subjectrange)
  betadbasis  <- create.bspline.basis(subjectrange,4)
  
  betadlist <- vector("list",pd)
  betadlist[[1]] <- betacbasis
  for (l in 2:pd) betadlist[[l]] <- betadbasis
  fRegressdList <- fRegress(subjecttempfd, xfddlist, betadlist)
  
  yhatdfdobj <- fRegressdList$yhatfdobj
  
  yhatdmat    <- eval.fd(1:10, yhatdfdobj$fd)
  ymatd       <- eval.fd(1:10, subjecttempfd)
  SSE <- sum((ymatd[,1:7] - yhatdmat[,1:7])^2)
  SSE
  
  ###########################################################################################################
  #######################################    Initialization   ###############################################
  ###########################################################################################################
  
  # hyperparameters
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
  
  sig2=0.00002
  sig0=0.00005
  sig2s<-c(sig2)
  sig0s<-c(sig0)
  
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
  BetaA<-c(beta)    
  BetaG<-c(beta_g)
  
  As<-c(a)
  Bs<-c(b)
  Cs<-c(c)
  thetas<-c(theta)
  
  s2g=1
  s2a=1
  S2A<-c(s2a)
  S2G<-c(s2g)
  # numCores <- detectCores()-1
  # numCores
  
  
  for (p in 2:M){
    #print(p)
    #source("/home/statgrads/krag57/Parallelized simulation5.R")
    XAbeta<-matrix(Xbind%*%beta,nrow = 7,byrow = T)
    XGbeta<-matrix(Xbind%*%beta_g,nrow = 7,byrow = T)
    
    BetaA = cbind(BetaA, beta)
    BetaG = cbind(BetaG, beta_g)
    
    pMu0<-muMu0(a,b,c,theta,d)
    
    registerDoParallel(10)
    Results=foreach(i=1:10) %dopar% {
      pos_DAlpGamMu0(y=y[2:3,i,],y0=y[1,i,],dalpha=resDAlpha[p-1,,i],dgamma=resDGamma[p-1,,i],
                     XAbeta=XAbeta[,i],XGbeta=XGbeta[,i],sigma=sqrt(sig2),betaASD=s2a,betaGSD=s2g)
    }
    #stopImplicitCluster()
    
    resDAlpha[p,,]<-do.call(cbind, lapply(Results, function(x) x[, 1]))
    resDGamma[p,,]<-do.call(cbind, lapply(Results, function(x) x[, 2]))
    
    resAlpha[p,1,]<-resDAlpha[p,1,]+1
    for (i in 2:7){
      resAlpha[p,i,]<-resAlpha[p,i-1,]+resDAlpha[p,i,]
    }
    
    resGamma[p,1,]<-resDGamma[p,1,]+1
    for (i in 2:7){
      resGamma[p,i,]<-resGamma[p,i-1,]+resDGamma[p,i,]
    }
    #print(p)
    sig2<-posSigma2(y,Alpha = t(resAlpha[p,,]),Gamma = t(resGamma[p,,]),Sigma2 = sig2,SSE=SSE)
    sig2s<-c(sig2s,sig2)
    
    #print(p)
    Mu0Sigma0<-posABCTheSig0(y0 = y[1,,],sig0 = sig0,a = a,b = b,c = c,theta = theta)
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
    
    capture.output(baa<-bridge(Xbind,matrix(log(resDAlpha[p,,]), ncol = 1),RJ=F,ab=c(2.92,cenAG*1.92)), file='NUL')
    beta<-colMeans(baa$beta[-(1:500),])
    s2a=(sqrt(mean(baa$s2)))
    
    capture.output(bag<-bridge(Xbind,matrix(log(resDGamma[p,,]), ncol = 1),RJ=F,ab=c(2.92,cenAG*1.92)), file='NUL')
    beta_g<-colMeans(bag$beta[-(1:500),])
    s2g=(sqrt(mean(bag$s2)))
    
    S2A<-c(S2A,s2a)
    S2G<-c(S2G,s2g)
  }
  
  AlphaHat<-array(0,c(7,5,dim(BetaA)[2]))
  for (r in 1:dim(BetaA)[2]){
    DA=matrix(0,nrow = 7,5)
    Sig0<-sqrt(sig0s[+r])
    DA<-matrix(exp(Xbind15[71:105,]%*%BetaA[,r]),nrow = 7)
    
    AThat=matrix(0,nrow = 7,5)
    AThat[1,]<-DA[1,]+1
    for (i in 2:7){
      AThat[i,]<-AThat[i-1,]+DA[i,]
    }
    AlphaHat[,,r]<-AThat
  }
  
  GammaHat<-array(0,c(7,5,dim(BetaG)[2]))
  for (r in 1:dim(BetaG)[2]){
    DG=matrix(0,nrow = 7,5)
    Sig0<-sqrt(sig0s[r])
    DG<-matrix(exp(Xbind15[71:105,]%*%BetaG[,r]),nrow = 7)
    
    GThat=matrix(0,nrow = 7,5)
    GThat[1,]<-DG[1,]+1
    for (j in 2:7){
      GThat[j,]<-GThat[j-1,]+DG[j,]
    }
    GammaHat[,,r]<-GThat
  }
  
  y0p<-Z[1,11:15,]
  r48<-Z[2,11:15,]
  r72<-Z[3,11:15,]
  
  y0pp<-matrix(as.vector(y0p),nrow = 7,ncol = 5,byrow = T)
  pre48<-array(0,c(5,7,dim(GammaHat)[3]))
  for (r in 1:dim(GammaHat)[3]){
    Sig<-(sig2s[r])
    p481dist<-rtmvnorm(1,mu =as.vector(y0pp*GammaHat[,,r]^(1-AlphaHat[,,r]^(-1))),sigma = diag(Sig,35), lb = rep(0,35),ub = rep(1,35))
    p481<-matrix(p481dist,5,7,byrow = T)
    pre48[,,r]<-p481
  }
  
  p48=apply(pre48, c(1,2), mean)
  
  pre72<-array(0,c(5,7,dim(GammaHat)[3]))
  for (r in 1:dim(GammaHat)[3]){
    Sig<-(sig2s[r])
    p721dist<-rtmvnorm(1,mu =as.vector(y0pp*GammaHat[,,r]^(1-AlphaHat[,,r]^(-2))),sigma = diag(Sig,35), lb = rep(0,35),ub = rep(1,35))
    p721<-matrix(p721dist,5,7,byrow = T)
    pre72[,,r]<-p721
  }
  
  p72=apply(pre72, c(1,2), mean)
  
  #################################################################################################
  #################################################################################################
  
  Al1=ParaInfoSubAlpha(AT[,1],resAlpha[(M/2):M,,],1)
  Al2=ParaInfoSubAlpha(AT[,2],resAlpha[(M/2):M,,],2)
  Al3=ParaInfoSubAlpha(AT[,3],resAlpha[(M/2):M,,],3)
  Al4=ParaInfoSubAlpha(AT[,4],resAlpha[(M/2):M,,],4)
  Al5=ParaInfoSubAlpha(AT[,5],resAlpha[(M/2):M,,],5)
  Al6=ParaInfoSubAlpha(AT[,6],resAlpha[(M/2):M,,],6)
  Al7=ParaInfoSubAlpha(AT[,7],resAlpha[(M/2):M,,],7)
  Al8=ParaInfoSubAlpha(AT[,8],resAlpha[(M/2):M,,],8)
  Al9=ParaInfoSubAlpha(AT[,9],resAlpha[(M/2):M,,],9)
  Al10=ParaInfoSubAlpha(AT[,10],resAlpha[(M/2):M,,],10)
  
  Ga1=ParaInfoSubGamma(GT[,1],resGamma,1)
  Ga2=ParaInfoSubGamma(GT[,2],resGamma,2)
  Ga3=ParaInfoSubGamma(GT[,3],resGamma,3)
  Ga4=ParaInfoSubGamma(GT[,4],resGamma,4)
  Ga5=ParaInfoSubGamma(GT[,5],resGamma,5)
  Ga6=ParaInfoSubGamma(GT[,6],resGamma,6)
  Ga7=ParaInfoSubGamma(GT[,7],resGamma,7)
  Ga8=ParaInfoSubGamma(GT[,8],resGamma,8)
  Ga9=ParaInfoSubGamma(GT[,9],resGamma,9)
  Ga10=ParaInfoSubGamma(GT[,10],resGamma,10)
  
  Si1=ParaInfoVectors(0.005,sqrt(sig0s),name="Sigma0")
  Si2=ParaInfoVectors(0.005,(sig2s),name="Sigma1")
  AA=ParaInfoVectors(at,(As),name="A")
  BB=ParaInfoVectors(bt,(Bs),name="B")
  CC=ParaInfoVectors(ct,(Cs),name="C")
  Th=ParaInfoVectors(thetat,(thetas),name="theta")
  SiB1=ParaInfoVectors(1,(S2A),name="PrVaAl")
  SiB2=ParaInfoVectors(0.5,(S2G),name="PrVaGa")
  
  BA=ParaInfoBetas(betaa,BetaA[,-(1:M/2)])
  BG=ParaInfoBetas(betag,BetaG[,-(1:M/2)])
  
  ComInfo<-rbind(Al1,Al2,Al3,Al4,Al5,Al6,Al7,Al8,Al9,Al10,
                 Ga1,Ga2,Ga3,Ga4,Ga5,Ga6,Ga7,Ga8,Ga9,Ga10,
                 Si1,Si2,AA,BB,CC,Th,SiB1,SiB2,BA,BG)
  
  
  Al1p<-ParaInfoSubAlphaHat(AT[,11],AlphaHat,1)
  Al2p<-ParaInfoSubAlphaHat(AT[,12],AlphaHat,2)
  Al3p<-ParaInfoSubAlphaHat(AT[,13],AlphaHat,3)
  Al4p<-ParaInfoSubAlphaHat(AT[,14],AlphaHat,4)
  Al5p<-ParaInfoSubAlphaHat(AT[,15],AlphaHat,5)
  
  Ga1p<-ParaInfoSubAlphaHat(GT[,11],GammaHat,1)
  Ga2p<-ParaInfoSubAlphaHat(GT[,12],GammaHat,2)
  Ga3p<-ParaInfoSubAlphaHat(GT[,13],GammaHat,3)
  Ga4p<-ParaInfoSubAlphaHat(GT[,14],GammaHat,4)
  Ga5p<-ParaInfoSubAlphaHat(GT[,15],GammaHat,5)
  
  
  y481p<-ParaInfoSubY(Z[2,11,],pre48,1)
  y482p<-ParaInfoSubY(Z[2,12,],pre48,2)
  y483p<-ParaInfoSubY(Z[2,13,],pre48,3)
  y484p<-ParaInfoSubY(Z[2,14,],pre48,4)
  y485p<-ParaInfoSubY(Z[2,15,],pre48,5)
  
  y721p<-ParaInfoSubY(Z[2,11,],pre72,1)
  y722p<-ParaInfoSubY(Z[2,12,],pre72,2)
  y723p<-ParaInfoSubY(Z[2,13,],pre72,3)
  y724p<-ParaInfoSubY(Z[2,14,],pre72,4)
  y725p<-ParaInfoSubY(Z[2,15,],pre72,5)
  
  PreInfo<-rbind(Al1p,Al2p,Al3p,Al4p,Al5p,
                 Ga1p,Ga2p,Ga3p,Ga4p,Ga5p,
                 y481p,y482p,y483p,y484p,y485p,
                 y721p,y722p,y723p,y724p,y725p
  )
  
  MAE48<-mean(abs(r48-p48))
  MAE72<-mean(abs(r72-p72))
  TMAE<-sum(MAE48,MAE72)/2
  
  PMSE48<-(mean((r48-p48)^2))
  PMSE72<-(mean((r72-p72)^2))
  TPMSE<-(sum(PMSE48,PMSE72)/2)
  
  PreEva<-c(MAE48,MAE72,TMAE,PMSE48,PMSE72,TPMSE)
  
  lisA=list("TAlpha"=AT,"TGamma"=GT,"TbetaA"=betaa,"TbetaG"=betag,"TrueHill"=ABCThetaTrue,
            "BetaA"=BetaA,"BetaG"=BetaG,"Alpha"=resAlpha,"Gamma"=resGamma,"Sigma0"=sqrt(sig0s),
            "Sigma1"=sig2s,"A0"=As,"B0"=Bs,"C0"=Cs,"Theta0"=thetas,"PrVaAl"=S2A,"PrVaGa"=S2G,
            "AlphaHat"=AlphaHat,"GammaHat"=GammaHat,"Y48"=pre48,"Y72"=pre72,"ParaInfo"=ComInfo,
            "PreInfo"=PreInfo,"Eval"=PreEva)
  
  save(lisA, file=paste0("/work/statgrads/krag57/Data",h,".RData"))
  
  EntireData <- append(EntireData, list(lisA))
  
  names(EntireData)[h] <- paste0("Data",h)
}
save(EntireData, file="/work/statgrads/krag57/EntireData10.RData")