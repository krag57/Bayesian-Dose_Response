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

# y = y[,i];y0 = y0[i];dalpha = dalpha[i];dgamma = dgamma[i];pre.dalpha = pre.dalpha[i];
# pre.dgamma = pre.dgamma[i];xabeta = XAbeta[i];sigma = sigma;betaSD = betaASD
priLikDAlpha<-function(y,y0,dalpha,dgamma,pre.dalpha,pre.dgamma,xabeta,sigma,betaSD){
  t<-1:2
  mu_tilde<-y0*(dgamma+pre.dgamma+1)^(1-(dalpha+pre.dalpha+1)^(-t))
  #dtnorm(y, low=0, upp=1, mean = mu_tilde, sd = sigma,log=T)
  likDAlpha<-sum(ifelse(dtnorm(y, low=0, upp=1, mean = mu_tilde, sd = sigma,log=T)==Inf,0,dtnorm(y, low=0, upp=1, mean = mu_tilde, sd = sigma,log=T)))
  #likDAlpha<-sum(dtnorm(y, low=0, upp=1, mean = mu_tilde, sd = sigma,log=T))#sum(log(dtruncnorm(y, a=0, b=1, mean = mu_tilde, sd = sigma)))dtnorm(y, low=0, upp=1, mean = mu_tilde, sd = sigma,log=T)
  priDAlpha<-log(dtruncnorm(dalpha, a=0,b=2, mean = xabeta, sd = betaSD))#sum(ifelse(dtnorm(dalpha, low=1, mean = xabeta, sd = betaSD,log=T)==Inf,0,dtnorm(dalpha, low=1, mean = xabeta, sd = betaSD,log=T)))#(dtnorm(dalpha, low=1, mean = xabeta, sd = betaSD,log=T))#dlnorm(dalpha,meanlog =xabeta,sdlog =betaSD,log=T)#dunif(dalpha,exp(xabeta)-(sqrt(3)*betaSD),max =exp(xabeta)+(sqrt(3)*betaSD),log=T)#
  return (likDAlpha+priDAlpha)
}

priLikDGamma<-function(y,y0,dalpha,dgamma,pre.dalpha,pre.dgamma,xgbeta,sigma,betaSD){
  t<-1:2
  mu_tilde<-y0*(dgamma+pre.dgamma+1)^(1-(dalpha+pre.dalpha+1)^(-t))
  likDGamma<-sum(ifelse(dtnorm(y, low=0, upp=1, mean = mu_tilde, sd = sigma,log=T)==Inf,0,dtnorm(y, low=0, upp=1, mean = mu_tilde, sd = sigma,log=T)))
  #likDGamma<-sum(dtnorm(y, low=0, upp=1, mean = mu_tilde, sd = sigma,log=T))#sum(log(dtruncnorm(y, a=0, b=1, mean = mu_tilde, sd = sigma)))
  priDGamma<-dtnorm(dgamma, low=0, upp=2, mean = xgbeta, sd = betaSD,log=T)#sum(ifelse(dtnorm(dgamma, low=1,  mean = xgbeta, sd = betaSD,log=T)==Inf,0,dtnorm(dgamma, low=1,  mean = xgbeta, sd = betaSD,log=T)))#dtnorm(dgamma, low=1,  mean = xgbeta, sd = betaSD,log=T)#log(dtruncnorm(dgamma, a=0, mean = xgbeta, sd = betaSD))#dlnorm(dgamma,meanlog =xgbeta,sdlog =betaSD,log=T)#dunif(dgamma,exp(xgbeta)-(sqrt(3)*betaSD),max =exp(xgbeta)+(sqrt(3)*betaSD),log=T)#
  return (likDGamma+priDGamma)
}

# y0=y[1,i,];y=y[2:3,i,];dalpha=resDAlpha[p-1,,i];dgamma=resDGamma[p-1,,i];
# XAbeta=XAbeta[,i];XGbeta=XGbeta[,i];sigma=sqrt(sig2);betaASD=s2a;betaGSD=s2g

pos_DAlpGamMu0<-function(y,y0,dalpha,dgamma,XAbeta,XGbeta,sigma,betaASD,betaGSD){
  pre.dalpha=c(0,dalpha)
  pre.dgamma=c(0,dgamma)
  for(i in 1:7){
    llik<-priLikDAlpha(y = y[,i],y0 = y0[i],dalpha = dalpha[i],dgamma = dgamma[i],pre.dalpha = pre.dalpha[i],pre.dgamma = pre.dgamma[i],xabeta = XAbeta[i],sigma = sigma,betaSD = betaASD)
    dalpha.s<-rtruncnorm(1,mean = dalpha[i],sd =.5,a=0,b=2)
    llikp<-priLikDAlpha(y = y[,i],y0 = y0[i],dalpha = dalpha.s,dgamma = dgamma[i],pre.dalpha = pre.dalpha[i],pre.dgamma = pre.dgamma[i],xabeta = XAbeta[i],sigma = sigma,betaSD = betaASD)
    r = exp(llikp - llik
            -log(dtruncnorm(dalpha.s, mean=dalpha[i], sd=.5, a=0,b=2))
            +log(dtruncnorm(dalpha[i], mean=dalpha.s, sd=.5, a=0,b=2)))
    accept<-min(1,r)
    z<-runif(1)
    if(z<accept){
      dalpha[i] = dalpha.s
      pre.dalpha[i+1]=dalpha.s
    }
    
    llikG<-priLikDGamma(y = y[,i],y0 = y0[i],dalpha=dalpha[i],dgamma = dgamma[i],pre.dalpha = pre.dalpha[i],pre.dgamma = pre.dgamma[i],xgbeta = XGbeta[i],sigma = sigma,betaSD = betaGSD)
    dgamma.s<-rtruncnorm(1,mean = dgamma[i],sd = .5, a=0,b=2)
    llikpG<-priLikDGamma(y = y[,i],y0 = y0[i],dalpha=dalpha[i],dgamma = dgamma.s,pre.dalpha = pre.dalpha[i],pre.dgamma = pre.dgamma[i],xgbeta = XGbeta[i],sigma = sigma,betaSD = betaGSD)
    r = exp(llikpG - llikG
            -log(dtruncnorm(dgamma.s, mean=dgamma[i], sd=.5, a=0,b=2))
            +log(dtruncnorm(dgamma[i], mean=dgamma.s, sd=.5, a=0,b=2)))
    accept<-min(1,r)
    z<-runif(1)
    if(z<accept){
      dgamma[i] = dgamma.s
      pre.dgamma[i+1]=dgamma.s
    }
  }
  return (cbind(dalpha,dgamma))
}


priLikSigma2<-function(y,Gamma,Alpha,Sigma2,SSE){
  likHoods<-array(0,c(2,10,7))
  y0<-y[1,,]
  for (k in 1:2){
    likHoods[k,,]<-y0*Gamma^(1-Alpha^(-k))
  }
  likY<-sum(ifelse(dtnorm(y[2:3,,], low=0, upp=1, mean = likHoods[1:2,,], sd = sqrt(Sigma2),log=T)==Inf,0,dtnorm(y[2:3,,], low=0, upp=1, mean = likHoods[1:2,,], sd = sqrt(Sigma2),log=T)))
  #likY<-sum(log(dtruncnorm(y[2:3,,], a=0, b=1, mean = likHoods[1:2,,], sd = sqrt(Sigma2))))#sum(msm::dtnorm(y[2:3,,],likHoods[1:2,,],sqrt(Sigma2),lower = 0,u=1,log = T))#
  priSigma2<-log(MCMCpack::dinvgamma(Sigma2,shape = 3,scale = SSE*2))#dinvgamma(10,shape = 3,scale = .019738e-05*2,log = T)
  return(likY+priSigma2)
}

#y=Y;Alpha = t(resAlpha[p,,]);Gamma = t(resGamma[p,,]);Sigma2 = sig2
posSigma2<-function(y,Gamma,Alpha,Sigma2,SSE){
  llik<-priLikSigma2(y = y,Gamma = Gamma,Alpha = Alpha,Sigma2 = Sigma2,SSE=SSE)
  Sigma2.s<-abs(rnorm(1,Sigma2,sd=0.0002))#rtruncnorm(1,mean =Sigma2,sd = 0.00005, a=0)#abs(rnorm(1,Sigma2,sd=0.005))#Sigma2.s=0.0002
  llikp<-priLikSigma2(y = y,Gamma = Gamma,Alpha = Alpha,Sigma2 = Sigma2.s,SSE=SSE)
  r = exp(llikp - llik)
  accept<-min(1,r)
  z<-runif(1)
  if(z<accept){
    Sigma2 = Sigma2.s
  }
  return (Sigma2)
}

muMu0<-function(a,b,c,theta,d) a+((b-a)/(1+(c/d)^theta))

priLikABCTheSig0<-function(y0,sig0,a,b,c,theta){
  d<-c(0.0032, 0.0100, 0.0316, 0.1000 ,0.3160 ,1.0000, 3.1600)
  mvmean<-muMu0(a,b,c,theta,d)
  mvsigma<-diag(sqrt(sig0),nrow = 7)
  likAll<-sum(dtmvnorm(y0, mu = mvmean,sigma = mvsigma,l=rep(0,7),u=rep(1,7),log = T))
  priSig0<-log(MCMCpack::dinvgamma(sig0,shape = 3,scale = 0.0001960908*2))#dgamma(sig0,shape = 6,scale = 0.008/5,log=T) # scale is determined from the main function.
  priA<-dunif(a,0,1,log = T)
  priB<-dunif(b,a,1,log=T)
  priC<-dunif(c,0,3.16,log=T)
  priTheta<-dunif(theta,0,3,log=T)
  return(likAll+priTheta+priSig0+priB+priA+priC)
}


posABCTheSig0<-function(y0,sig0,a,b,c,theta){
  likeH<-priLikABCTheSig0(y0 = y0,sig0 = sig0,a = a,b = b,c = c,theta = theta)
  a.s<-rtruncnorm(1,a,sd=0.05,a = 0,b=1)
  likeHS<-priLikABCTheSig0(y0 = y0,sig0 = sig0,a = a.s,b = b,c = c,theta = theta)
  r = exp(likeHS - likeH
          -log(dtruncnorm(a.s,mean =a,sd = 0.05, a=0,b=1))
          +log(dtruncnorm(a,mean =a.s,sd = 0.05, a=0,b=1)))
  accept<-min(1,r)
  z<-runif(1)
  if(z<accept){
    a=a.s
  }
  
  likeH<-priLikABCTheSig0(y0 = y0,sig0 = sig0,a = a,b = b,c = c,theta = theta)
  b.s<-rtruncnorm(1,mean = b,sd=0.05,a = a,b=1)
  likeHS<-priLikABCTheSig0(y0 = y0,sig0 = sig0,a = a,b = b.s,c = c,theta = theta)
  r = exp(likeHS - likeH
          -log(dtruncnorm(x = b.s,mean =b,sd = 0.05, a=a, b=1))
          +log(dtruncnorm(x = b,mean =b.s,sd = 0.05, a=a,  b=1)))
  accept<-min(1,r)
  z<-runif(1)
  if(z<accept){
    b=b.s
  }
  
  likeH<-priLikABCTheSig0(y0 = y0,sig0 = sig0,a = a,b = b,c = c,theta = theta)
  c.s<-rtruncnorm(1,mean=c,sd=0.5,a = 0,b = 3.16)
  likeHS<-priLikABCTheSig0(y0 = y0,sig0 = sig0,a = a,b = b,c = c.s,theta = theta)
  r = exp(likeHS - likeH
          -log(dtruncnorm(x = c.s,mean =c,sd = 0.5, a=0,b = 3.16))
          +log(dtruncnorm(x = c,mean =c.s,sd = 0.5, a=0, b = 3.16)))
  accept<-min(1,r)
  z<-runif(1)
  if(z<accept){
    c=c.s
  }
  
  likeH<-priLikABCTheSig0(y0 = y0,sig0 = sig0,a = a,b = b,c = c,theta = theta)
  theta.s<-rtruncnorm(1,mean=theta,sd=0.5,a = 0,b=3)
  likeHS<-priLikABCTheSig0(y0 = y0,sig0 = sig0,a = a,b = b,c = c,theta = theta.s)
  r = exp(likeHS - likeH
          -log(dtruncnorm(theta.s,mean =theta,sd = 0.5, a=0, b = 3))
          +log(dtruncnorm(theta,mean =theta.s,sd = 0.5, a=0, b=3)))
  accept<-min(1,r)
  z<-runif(1)
  if(z<accept){
    theta=theta.s
  }
  
  likeH<-priLikABCTheSig0(y0 = y0,sig0 = sig0,a = a,b = b,c = c,theta = theta)
  sig0.s<-abs(rnorm(1,sig0,sd=0.00005))#rtruncnorm(1,mean =sig0,sd = 0.005, a=0)#
  likeHS<-priLikABCTheSig0(y0 = y0,sig0 = sig0.s,a = a,b = b,c = c,theta = theta)
  r = exp(likeHS - likeH)
  #-log(dtruncnorm(sig0.s,mean =sig0,sd = 0.005, a=0))
  #+log(dtruncnorm(sig0,mean =sig0.s,sd = 0.005, a=0)))
  accept<-min(1,r)
  z<-runif(1)
  if(z<accept){
    sig0=sig0.s
  }
  
  return (c(a,b,c,theta,sig0))
}

##########################################################################
############################# PLots ######################################
##########################################################################

alphaGamma<-function(arr,dose,subject){
  count=c()
  for(i in 1:dim(arr)[1]){
    count=c(count,arr[i,dose,subject])
  }
  return (count)
}

alphaGammaHat<-function(arr,dose,subject){
  count=c()
  for(i in 1:dim(arr)[3]){
    count=c(count,arr[dose,subject,i])
  }
  return (count)
}

yHat<-function(arr,dose,subject){
  count=c()
  for(i in 1:dim(arr)[3]){
    count=c(count,arr[subject,dose,i])
  }
  return (count)
}

ParaInfoSubY<-function(tv,arr,subject){
  #par(mfrow=c(2,5))
  info<-c()
  for (i in 1:7){
    pr=round(mean(yHat(arr,i,subject)),3)
    prl=round(quantile(yHat(arr,i,subject),0.025),3)
    pru=round(quantile(yHat(arr,i,subject),0.975),3)
    #ge=round(geweke.diag(as.mcmc(alphaGammaHat(arr,dose,i)),frac2 = .1)$z,3)
    #info=rbind(info,c(round(tv[i],3),prl,pr,pru,ge))
    info=rbind(info,c(round(tv[i],3),prl,pr,pru))
  }
  prmatrix(info,rowlab = rep("",7))
  #return(info)
  #write.csv(info,file = paste0("~/Documents/GitHub/Bayesian-Dose_Response/Alpha",dose,".csv"))
}

ParaInfoSubAlphaHat<-function(tv,arr,dose){
  #par(mfrow=c(2,5))
  info<-c()
  for (i in 1:7){
    pr=round(mean(alphaGammaHat(arr,i,dose)),3)
    prl=round(quantile(alphaGammaHat(arr,i,dose),0.025),3)
    pru=round(quantile(alphaGammaHat(arr,i,dose),0.975),3)
    #ge=round(geweke.diag(as.mcmc(alphaGammaHat(arr,dose,i)),frac2 = .1)$z,3)
    #info=rbind(info,c(round(tv[i],3),prl,pr,pru,ge))
    info=rbind(info,c(round(tv[i],3),prl,pr,pru))
  }
  prmatrix(info,rowlab = rep("",7))
  #return(info)
  #write.csv(info,file = paste0("~/Documents/GitHub/Bayesian-Dose_Response/Alpha",dose,".csv"))
}





tracePlotSub<-function(arr,dose){
  par(mfrow=c(2,5))
  for (i in 1:10){
    traceplot(as.mcmc(alphaGamma(arr,dose,i)))
  }
}


ParaInfoSubAlpha<-function(tv,arr,dose){
  #par(mfrow=c(2,5))
  info<-c()
  for (i in 1:7){
    pr=round(mean(alphaGamma(arr,i,dose)),3)
    prl=round(quantile(alphaGamma(arr,i,dose),0.025),3)
    pru=round(quantile(alphaGamma(arr,i,dose),0.975),3)
    ge=round(geweke.diag(as.mcmc(alphaGamma(arr,i,dose)),frac2 = .1)$z,3)
    info=rbind(info,c(round(tv[i],3),prl,pr,pru,ge))
  }
  prmatrix(info,rowlab = rep("",7))
  #return(info)
  #write.csv(info,file = paste0("~/Documents/GitHub/Bayesian-Dose_Response/Alpha",dose,".csv"))
}

ParaInfoSubGamma<-function(tv,arr,dose){
  #par(mfrow=c(2,5))
  info<-c()
  for (i in 1:7){
    pr=round(mean(alphaGamma(arr,dose,i)),3)
    prl=round(quantile(alphaGamma(arr,dose,i),0.025),3)
    pru=round(quantile(alphaGamma(arr,dose,i),0.975),3)
    ge=round(geweke.diag(as.mcmc(alphaGamma(arr,dose,i)),frac2 = .1)$z,3)
    info=rbind(info,c(round(tv[i],3),prl,pr,pru,ge))
  }
  prmatrix(info,rowlab = rep("",7))
  #return(info)
  #write.csv(info,file = paste0("~/Documents/GitHub/Bayesian-Dose_Response/Gamma",dose,".csv"))
}

ParaInfoVectors<-function(tv,arr,name){
  info<-c()
  pr=round(mean(arr),3)
  prl=round(quantile(arr,0.025),3)
  pru=round(quantile(arr,0.975),3)
  ge=round(geweke.diag(as.mcmc(arr),frac2 = .1)$z,3)
  info=rbind(info,c(round(tv,3),prl,pr,pru,ge))
  return(info)
  #write.csv(info,file = paste0("~/Documents/GitHub/Bayesian-Dose_Response/",name,".csv"))
}

ParaInfoBetas<-function(be,arr){
  sumBetaA<-cbind(be,
                  round(apply(arr,1,function (x) quantile(x,0.025)),3),
                  round(rowMeans(arr),3),
                  round(apply(arr,1,function (x) quantile(x,0.975)),3),
                  round(geweke.diag(as.mcmc(t(arr)))$z,3))
  
  #write.csv(sumBetaA,file = paste0("/work/statgrads/krag57/AZ",j,"/Summary beta alpha.csv"))
  prmatrix(sumBetaA,rowlab = rep("",14))
  #return(sumBetaA)
  #write.csv(info,file = paste0("~/Documents/GitHub/Bayesian-Dose_Response/",name,".csv"))
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

gplotrow<-function(data,tru,name="Alpha 1"){
  nrow=dim(data)[1]
  x=1:nrow
  low=apply(data,1,function (x) quantile(x,0.025))
  upp=apply(data,1,function (x) quantile(x,0.975))
  pre=rowMeans(data)
  
  pl=data.frame(x=x,y=pre,z=tru,
                l=low,
                u=upp)
  
  ggplot(pl) + 
    geom_line(aes(x,y,colour="blue")) +
    geom_point(aes(x,z,colour="darkred")) +
    theme_gray() +
    labs(title = paste0("Prediction including 95%CI ",name),x="Dose",y="Values") +
    geom_ribbon(aes(x,ymin=l, ymax=u,color="lightblue"),fill="lightblue", alpha=0.3)+
    scale_x_continuous(breaks=x)+
    scale_color_manual(name=NULL,values = c("blue", "red","lightblue"),labels=c("Predicted", "True","95% CI"))+
    theme(plot.title = element_text(hjust=0.5),legend.position = c(0.95, 0.95),legend.justification = c("right", "top"),legend.background = element_rect(fill=NA))
  
  #ggsave(paste0(name,".png"))
}

gplotrow1<-function(data,tru,name="Alpha 1"){
  nrow=dim(data)[1]
  x=1:nrow
  low=apply(data,1,function (x) quantile(x,0.025))
  upp=apply(data,1,function (x) quantile(x,0.975))
  pre=rowMeans(data)
  
  pl=data.frame(x=x,y=pre,z=tru,
                l=low,
                u=upp)
  
  ggplot(pl) + 
    geom_line(aes(x,y,colour="blue")) +
    geom_point(aes(x,z,colour="darkred")) +
    theme_gray() +
    labs(title = paste0("Prediction including 95%CI ",name),x="Predictors",y="Values") +
    geom_ribbon(aes(x,ymin=l, ymax=u,color="lightblue"),fill="lightblue", alpha=0.3)+
    scale_x_continuous(breaks=x)+
    scale_y_continuous(limits=c(-6,6))+
    scale_color_manual(name=NULL,values = c("blue", "red","lightblue"),labels=c("Predicted", "True","95% CI"))+
    theme(plot.title = element_text(hjust=0.5),legend.position = c(0.95, 0.95),legend.justification = c("right", "top"),legend.background = element_rect(fill=NA))
  
  #ggsave(paste0(name,".png"))
}

gplotpred1<-function(pre,low,upp,tru,name="Alpha 1",posi="none"){
  nrow=length(pre)
  x=1:nrow
  pl=data.frame(x=x,y=pre,z=tru,l=low,u=upp)
  ggplot(pl) + 
    geom_line(aes(x,y,colour="blue")) +
    geom_point(aes(x,z,colour="darkred")) +
    theme_gray() +
    labs(title = name,x="Dose",y="Values") +
    geom_ribbon(aes(x,ymin=l, ymax=u,color="lightblue"),fill="lightblue", alpha=0.3)+
    scale_x_continuous(breaks=x)+
    scale_color_manual(name=NULL,values = c("blue", "red","lightblue"),labels=c("Predicted", "True","95% CI"))+
    theme(legend.position = posi,legend.justification = c("right", "top"),legend.background = element_rect(fill=NA),
          plot.title = element_text(hjust=0.5,size=10),axis.title.x = element_text(size=8),
          axis.title.y = element_text(size=8),axis.text.x = element_text(size=7),
          axis.text.y = element_text(size=7),legend.text=element_text(size=8),legend.key.size = unit(.75, "cm"),
          legend.key.width = unit(0.75,"cm") )
  
  #ggsave(paste0(name,".png"))
}

gplotpred<-function(df,dfl,dfu,dft,variable="AZ-638 at t=48"){
  n=dim(df)[1]
  plots <- vector("list", n)
  for (i in 1:n){
    plots[[i]] <- gplotpred1(df[i,],dfl[i,],dfu[i,],dft[i,],name=paste0("Subject ",i))
    p2=get_legend(gplotpred1(df[i,],dfl[i,],dfu[i,],dft[i,],name=paste0("Subject ",i),posi=c(0.75, 0.75)))
  }
  plots[[i+1]]=p2
  
  figure=do.call(ggarrange, plots)#c(plots,common.legend = TRUE, legend="right"))
  # theme(plot.subtitle = element_text(hjust=0.5,size=8),axis.title.x = element_text(size=8),
  #       axis.title.y = element_text(size=8),axis.text.x = element_text(size=5),
  #       axis.text.y = element_text(size=5))
  annotate_figure(figure,top = text_grob(paste0("Prediction of  ",variable),face = "bold",size = 16))
  #ggsave(paste0("Trace plot of ",variable,".png"))
}

gplotpred2<-function(pre,low,upp,tru,name="Alpha 1",posi="none"){
  nrow=length(pre)
  x=1:nrow
  pl=data.frame(x=x,y=pre,z=tru,l=low,u=upp)
  ggplot(pl) + 
    geom_line(aes(x,y,colour="blue")) +
    geom_point(aes(x,z,colour="darkred")) +
    theme_gray() +
    labs(title = name,x="Doses",y="Values") +
    geom_ribbon(aes(x,ymin=l, ymax=u,color="lightblue"),fill="lightblue", alpha=0.3)+
    scale_x_continuous(breaks=x)+
    scale_color_manual(name=NULL,values = c("blue", "red","lightblue"),labels=c("Predicted", "True","95% CI"))+
    theme(legend.position = posi,legend.justification = c("right", "top"),legend.background = element_rect(fill=NA),
          plot.title = element_text(hjust=0.5,size=10),axis.title.x = element_text(size=8),
          axis.title.y = element_text(size=8),axis.text.x = element_text(size=7),
          axis.text.y = element_text(size=7),legend.text=element_text(size=8),legend.key.size = unit(.75, "cm"),
          legend.key.width = unit(0.75,"cm") )
  
  #ggsave(paste0(name,".png"))
}

gplotpred3<-function(df,dfl,dfu,dft,variable="Alpha"){
  n=dim(df)[1]
  plots <- vector("list", n)
  for (i in 1:n){
    plots[[i]] <- gplotpred2(df[i,],dfl[i,],dfu[i,],dft[i,],name=paste0("Subject ",i))
    p2=get_legend(gplotpred2(df[i,],dfl[i,],dfu[i,],dft[i,],name=paste0("Subject ",i),posi=c(0.75, 0.75)))
  }
  plots[[i+1]]=p2
  
  figure=do.call(ggarrange, plots)#c(plots,common.legend = TRUE, legend="right"))
  # theme(plot.subtitle = element_text(hjust=0.5,size=8),axis.title.x = element_text(size=8),
  #       axis.title.y = element_text(size=8),axis.text.x = element_text(size=5),
  #       axis.text.y = element_text(size=5))
  annotate_figure(figure,top = text_grob(paste0("Prediction including 95%CI ",variable,"s"),face = "bold",size = 16))
  #ggsave(paste0("Trace plot of ",variable,".png"))
}

gplotcol<-function(data,tru,name="Alpha 1"){
  
  ncol=dim(data)[2]
  x=1:ncol
  low=apply(data,2,function (x) quantile(x,0.025))
  upp=apply(data,2,function (x) quantile(x,0.975))
  pre=colMeans(data)
  
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
  
  #ggsave(paste0(name,".png"))
}

# gplot(1:15,rowMeans(temp_a1[,-(1:M/2)]),apply(temp_a1[,-(1:M/2)],1,function (x) quantile(x,0.025)),
#       apply(temp_a1[,-(1:M/2)],1,function (x) quantile(x,0.975)),AT[,1])
# 
# df=temp_a1[,-(1:M/2)]
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
  annotate_figure(figure,top = text_grob(paste0("Trace plot of ",variable),face = "bold"))
  #ggsave(paste0("Trace plot of ",variable,".png"))
}
#tracegg(df)
tracesingle<-function(df,variable="Sigma0"){
  m=length(df)
  ggplot(data.frame(df),aes(x=1:m,y=df)) +
    geom_line()+
    labs(title = paste0("Trace plot of ",variable),x="Iteration",y="Values")+
    #scale_x_continuous(breaks=c(1,m/2,m))+
    #scale_y_continuous(breaks=as.character(round(min(df[i,])-1,2)),as.character(round(max(df[i,])+1,2)),as.character(round(mean(df[i,]),2)))+
    theme(plot.title = element_text(hjust=0.5,size=16,face="bold"),axis.title.x = element_text(size=10),
          axis.title.y = element_text(size=10),axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))
}


