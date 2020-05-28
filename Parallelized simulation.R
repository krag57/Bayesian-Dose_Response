integral<-function(cova=t(can_c32_AZ_1)[,2]){
  tempbasis<-create.bspline.basis(c(1,5),nbasis = 3,norder = 3)
  tempsmooth1<-smooth.basis(seq(1,5,le=3),cova,tempbasis)
  xfdj       <- tempsmooth1$fd
  xbasis     <- xfdj$basis
  xnbasis    <- xbasis$nbasis
  xrng       <- xbasis$rangeval
  nfine      <- 5000
  tfine      <- seq(xrng[1], xrng[2], len=nfine)
  deltat     <- tfine[2]-tfine[1]
  xmat       <- eval.fd(tfine, xfdj)
  fitj       <- deltat*(crossprod(xmat,rep(1,nfine)) -0.5*(outer(xmat[1,],1) +outer(xmat[nfine,],1)))
  return (abs(fitj))
}

designCell<-function(data=AZ_628,cell="C32"){
  x10<-matrix(0,14,5)
  colnames(x10)=as.character(unique(data$Time.Point..hr.))
  sub<-unique(data$Drug.Concentration..uM.)
  x.0<-list()
  X0<-matrix(0,7,14)
  for (j in 1:7){
    for (i in 1:14){
      x10[i,]<-data[(data$Drug.Concentration..uM. ==sub[j] & data$Cell.Line.Name ==cell),(i+4)]
    }
    x.0[[j]]<-x10[,3:5]
    X0[j,]<-as.vector(integral(t(x.0[[j]])))
  }
  return (X0)
}

priLikDAlpha<-function(y,y0,dalpha,dgamma,pre.dalpha,xabeta,sigma,betaSD){
  t<-1:2
  mu_tilde<-y0*(dgamma+pre.dalpha+1)^(1-(dalpha+pre.dalpha+1)^(-t))
  likDAlpha<-sum(log(dtruncnorm(y, a=0, b=1, mean = mu_tilde, sd = sigma)))
  priDAlpha<-dlnorm(dalpha,meanlog =xabeta,sdlog =betaSD,log=T)
  return (likDAlpha+priDAlpha)
}

priLikDGamma<-function(y,y0,dalpha,dgamma,pre.dgamma,xgbeta,sigma,betaSD){
  t<-1:2
  mu_tilde<-y0*(dgamma+pre.dgamma+1)^(1-(dalpha+pre.dgamma+1)^(-t))
  likDGamma<-sum(log(dtruncnorm(y, a=0, b=1, mean = mu_tilde, sd = sigma)))
  priDGamma<-dlnorm(dgamma,meanlog =xgbeta,sdlog =betaSD,log=T)
  return (likDGamma+priDGamma)
}

priLikDMu0<-function(y0,mu0,sigma,pMu0,pSigma0){
  likMu0<-log(dtruncnorm(y0, a=0, b=1, mean = mu0, sd = (sigma)))
  priMu0<-log(dtruncnorm(mu0,mean = pMu0,sd = pSigma0,a = 0,b=1))
  return (likMu0+priMu0)
}

pos_DAlpGamMu0<-function(y,y0,mu0,dalpha,dgamma,XAbeta,XGbeta,sigma,betaASD,betaGSD,pMu0,pSigma0){
  pre.dalpha=c(0,dalpha)
  pre.dgamma=c(0,dgamma)
  for(i in 1:7){
    llik<-priLikDAlpha(y = y[,i],y0 = y0[i],dalpha = dalpha[i],dgamma = dgamma[i],pre.dalpha = pre.dalpha[i],xabeta = XAbeta[i],sigma = sigma,betaSD = betaASD)
    dalpha.s<-rtruncnorm(1,mean = dalpha[i],sd = .1,a=0)
    llikp<-priLikDAlpha(y = y[,i],y0 = y0[i],dalpha = dalpha.s,dgamma = dgamma[i],pre.dalpha = pre.dalpha[i],xabeta = XAbeta[i],sigma = sigma,betaSD = betaASD)
    r = exp(llikp - llik
            -log(dtruncnorm(dalpha.s, mean=dalpha[i], sd=.1, a=0))
            +log(dtruncnorm(dalpha[i], mean=dalpha.s, sd=.1, a=0)))
    accept<-min(1,r)
    z<-runif(1)
    if(z<accept){
      dalpha[i] = dalpha.s
      pre.dalpha[i+1]=dalpha.s
    }
    
    llikG<-priLikDGamma(y = y[,i],y0 = y0[i],dalpha=dalpha[i],dgamma = dgamma[i],pre.dgamma = pre.dgamma[i],xgbeta = XGbeta[i],sigma = sigma,betaSD = betaGSD)
    dgamma.s<-rtruncnorm(1,mean = dgamma[i],sd = .1, a=0)
    llikpG<-priLikDGamma(y = y[,i],y0 = y0[i],dalpha=dalpha[i],dgamma = dgamma.s,pre.dgamma = pre.dgamma[i],xgbeta = XGbeta[i],sigma = sigma,betaSD = betaGSD)
    r = exp(llikpG - llikG
            -log(dtruncnorm(dgamma.s, mean=dgamma[i], sd=.1, a=0))
            +log(dtruncnorm(dgamma[i], mean=dgamma.s, sd=.1, a=0)))
    accept<-min(1,r)
    z<-runif(1)
    if(z<accept){
      dgamma[i] = dgamma.s
      pre.dgamma[i+1]=dgamma.s
    }
    
    llikM<-priLikDMu0(y0 = y0[i],mu0 = mu0[i],sigma = sigma,pMu0 = pMu0[i],pSigma0 = pSigma0)
    Mu0.s<-rtruncnorm(1,mean =mu0[i],sd = 0.05, a=0,b=1)
    llikpM<-priLikDMu0(y0 = y0[i],mu0 = Mu0.s,sigma = sigma,pMu0 = pMu0[i],pSigma0 = pSigma0)
    r = exp(llikpM - llikM
            -log(dtruncnorm(Mu0.s,mean =mu0[i],sd = 0.05, a=0, b=1))
            +log(dtruncnorm(mu0[i],mean =Mu0.s,sd = 0.05, a=0, b=1)))
    accept<-min(1,r)
    z<-runif(1)
    if(z<accept){
      mu0[i] = Mu0.s
    }
  }
  return (cbind(dalpha,dgamma,mu0))
}

alphaGamma<-function(arr,dose,subject){
  count=c()
  for(i in 1:dim(arr)[1]){
    count=c(count,arr[i,dose,subject])
  }
  return (count)
}

trplot<-function(x){
  par(mfrow=c(2,5))
  for(l in 1:10){
    traceplot(as.mcmc(x[l,-c(1:1000)]))
  }
  title(deparse(substitute(x)), line = -2, outer = TRUE)
}

colMode<-function(x){apply(x,2,estimate_mode)}
rowMode<-function(x){apply(x,1,estimate_mode)}

estimate_mode<-function(x){
  d <- density(x)
  return (round(d$x[which.max(d$y)],6))
}

priLikSigma2<-function(y,mu0,Gamma,Alpha,Sigma2){
  likHoods<-array(0,c(2,10,7))
  y0<-y[1,,]
  for (k in 1:2){
    likHoods[k,,]<-y0*Gamma^(1-Alpha^(-k))
  }
  likY<-sum(log(dtruncnorm(y[2:3,,], a=0, b=1, mean = likHoods, sd = sqrt(Sigma2))))
  likY0<-sum(log(dtruncnorm(y0, a=0, b=1, mean = mu0, sd = sqrt(Sigma2))))
  priSigma2<-dinvgamma(Sigma2,shape = 2,scale = 1,log = T)
  return(likY+likY0+priSigma2)
}

posSigma2<-function(y,mu0,Gamma,Alpha,Sigma2){
  llik<-priLikSigma2(y = y,mu0 = mu0,Gamma = Gamma,Alpha = Alpha,Sigma2 = Sigma2)
  Sigma2.s<-abs(rnorm(1,Sigma2,sd=0.005))
  llikp<-priLikSigma2(y = y,mu0 = mu0,Gamma = Gamma,Alpha = Alpha,Sigma2 = Sigma2.s)
  r = exp(llikp - llik)
  accept<-min(1,r)
  z<-runif(1)
  if(z<accept){
    Sigma2 = Sigma2.s
  }
  return (Sigma2)
}


muMu0<-function(a,b,c,theta,d) a+((b-a)/(1+(c/d)^theta))
priLikABCTheSig0<-function(mu0,sig0,a,b,c,theta){
  d<-c(0.0032, 0.0100, 0.0316, 0.1000 ,0.3160 ,1.0000, 3.1600)
  mvmean<-muMu0(a,b,c,theta,d)
  mvsigma<-diag(sqrt(sig0),nrow = 7)
  likAll<-sum(dtmvnorm(mu0, mu = mvmean,sigma = mvsigma,l=rep(0,7),u=rep(1,7),log = T))
  priSig0<-dinvgamma(sig0,shape=1,scale=1,log=T)
  priA<-dunif(a,0,1,log = T)
  priB<-dunif(b,a,1,log=T)
  priC<-dunif(c,0,3.16,log=T)
  priTheta<-dunif(theta,0,3,log=T)
  return(likAll+priTheta+priSig0+priB+priA+priC)
}


posABCTheSig0<-function(mu0,sig0,a,b,c,theta){
  likeH<-priLikABCTheSig0(mu0 = mu0,sig0 = sig0,a = a,b = b,c = c,theta = theta)
  a.s<-rtruncnorm(1,a,sd=0.05,a = 0,b=1)
  likeHS<-priLikABCTheSig0(mu0 = mu0,sig0 = sig0,a = a.s,b = b,c = c,theta = theta)
  r = exp(likeHS - likeH
          -log(dtruncnorm(a.s,mean =a,sd = 0.05, a=0,b=1))
          +log(dtruncnorm(a,mean =a.s,sd = 0.05, a=0,b=1)))
  accept<-min(1,r)
  z<-runif(1)
  if(z<accept){
    a=a.s
  }
  
  likeH<-priLikABCTheSig0(mu0 = mu0,sig0 = sig0,a = a,b = b,c = c,theta = theta)
  b.s<-rtruncnorm(1,mean = b,sd=0.05,a = a,b=1)
  likeHS<-priLikABCTheSig0(mu0 = mu0,sig0 = sig0,a = a,b = b.s,c = c,theta = theta)
  r = exp(likeHS - likeH
          -log(dtruncnorm(x = b.s,mean =b,sd = 0.05, a=a, b=1))
          +log(dtruncnorm(x = b,mean =b.s,sd = 0.05, a=a,  b=1)))
  accept<-min(1,r)
  z<-runif(1)
  if(z<accept){
    b=b.s
  }
  
  likeH<-priLikABCTheSig0(mu0 = mu0,sig0 = sig0,a = a,b = b,c = c,theta = theta)
  c.s<-rtruncnorm(1,mean=c,sd=0.5,a = 0,b = 3.16)
  likeHS<-priLikABCTheSig0(mu0 = mu0,sig0 = sig0,a = a,b = b,c = c.s,theta = theta)
  r = exp(likeHS - likeH
          -log(dtruncnorm(x = c.s,mean =c,sd = 0.5, a=0,b = 3.16))
          +log(dtruncnorm(x = c,mean =c.s,sd = 0.5, a=0, b = 3.16)))
  accept<-min(1,r)
  z<-runif(1)
  if(z<accept){
    c=c.s
  }
  
  likeH<-priLikABCTheSig0(mu0 = mu0,sig0 = sig0,a = a,b = b,c = c,theta = theta)
  theta.s<-rtruncnorm(1,mean=theta,sd=0.5,a = 0,b=3)
  likeHS<-priLikABCTheSig0(mu0 = mu0,sig0 = sig0,a = a,b = b,c = c,theta = theta.s)
  r = exp(likeHS - likeH
          -log(dtruncnorm(theta.s,mean =theta,sd = 0.5, a=0, b = 3))
          +log(dtruncnorm(theta,mean =theta.s,sd = 0.5, a=0, b=3)))
  accept<-min(1,r)
  z<-runif(1)
  if(z<accept){
    theta=theta.s
  }
  
  likeH<-priLikABCTheSig0(mu0 = mu0,sig0 = sig0,a = a,b = b,c = c,theta = theta)
  sig0.s<-abs(rnorm(1,sig0,sd=0.05))
  likeHS<-priLikABCTheSig0(mu0 = mu0,sig0 = sig0.s,a = a,b = b,c = c,theta = theta)
  r = exp(likeHS - likeH)
  accept<-min(1,r)
  z<-runif(1)
  if(z<accept){
    sig0=sig0.s
  }
  
  return (c(a,b,c,theta,sig0))
}




