#sampling function
rsample_ns=function(n,nt,nv,ntest,Betay,random,c,cbeta){
  pmean=rep(5,3)
  psigma=10*diag(3)
  X1=mvrnorm(n, pmean, psigma)
  X2=cbind(rbinom(n,1,0.5),rbinom(n,1,0.5),rbinom(n,1,0.5))
  X3=cbind(logit(runif(n)),logit(runif(n)),logit(runif(n)))
  x0=rep(1,n)
  nx=as.matrix(cbind(x0,X1,X2,X3))
  p=plogis(nx%*%Betay)
  #simulate y
  y=rbinom(n,1,p)
  #prevalence of y in the population  
  prevalence=mean(p) 
  #simulate s
  s=rep(0,n)
  if (random==TRUE) {
    for (i in 1:n) {
      if (y[i]==1){
        s[i]=rbinom(1,1,c)} else {s[i]==0}
    }
  }
  if (random==FALSE) {
    c1=plogis(nx%*%cbeta)
    for (i in 1:n) {
      if (y[i]==1){
        s[i]=rbinom(1,1,c1[i])} else {s[i]==0}
    }
  }
  pop=data.frame(nx,y,s,p)   ##population
  
  train_int=sample.int(nrow(pop), size=nt , replace = FALSE, prob = NULL)
  train=pop[train_int,]
  
  ###testing sample
  pop1=pop[-train_int,]
  test_int=sample.int(nrow(pop1), size=ntest , replace = FALSE, prob = NULL)
  test=pop1[test_int,]# sample ntest as testing sample
  
  ###validation sample from whole population
  pop2=pop1[-test_int,]
  valid_int=sample.int(nrow(pop2), size=nv , replace = FALSE, prob = NULL)
  valid=pop2[valid_int,] #validation sample of whole population

  return(list(prevalence=prevalence,train=train,valid=valid,test=test))
}

####optim function

LZ_c_optim <- function(train,iter, tol,Betay){
  
  ##likelihood function
  tx=as.matrix(train[,1:length(Betay)])
  nn=length(Betay)
  
  lik <- function(par) {
    beta <- par[1:nn]
    c<- par[nn+1]
    #lamda1 <- par[nn+2]
    
    sum1 <- -sum(train$s*log(plogis(tx %*% beta)*plogis(c))+(1-train$s)*log(1-plogis(c)*plogis(tx %*% beta)))
    #sum2 <- -sum(tx[,1]*log(1/lamda1))
    
    sum1
  }

  ctrl <- list(outer.iter=iter,tol=tol,trace=0)
  
  ##starting value
  m <- glm(train$s ~ tx[,-1], family=binomial())
  betahat <- m$coefficients
  chat <- sum(train$s)/length(train$s)
  lchat=log(chat/(1-chat))
  par0 <- c(betahat, lchat)
  
  ### optim optimization
  fit <-optim(par0, fn=lik, hessian = T,method = "BFGS") 
  theta=fit$par
  theta[length(par0)]=plogis(fit$par[length(par0)])
  if (fit$convergence==0) {conv=1} else {conv=0} # conv=1 if converged
  hess <- fit$hessian
  cov=solve(hess)
  se=sqrt(diag(cov))
  cc=theta[length(par0)]
  delta=cc*(1-cc)
  se[length(par0)]=se[length(par0)]*delta
  bcov=cov[1:length(Betay),1:length(Betay)]
  
  #delta method to compute ase of q
  betay=theta[1:length(Betay)]
  aa=as.vector(plogis(tx %*% betay)*(1-plogis(tx %*% betay)))
  gg=matrix(0,nrow=nrow(tx),ncol=ncol(tx))
  for (i in 1:nrow(tx)){
    gg[i,]=aa[i]*tx[i,]
  }
  gbeta=apply(gg,2,mean)
  gbeta=as.matrix(gbeta,nrow=length(gbeta))
  qcov_b=t(gbeta)%*%bcov%*%gbeta
  
  #empirical on F
  fvec=plogis(tx %*% betay)-mean(plogis(tx %*% betay))
  qcov_f=var(fvec)/length(fvec)
  
  #SE of qhat
  qcov=qcov_b+qcov_f
  qse=sqrt(qcov)
  
  return(list(theta=theta,se=se,conv=conv,cov=cov,qse=qse,bcov=bcov))
}



nonpara_ab<-function(s, predy, v, pp,cc1){
  yhat <- predy
  value <- numeric(length(v))
  for (j in 1:(length(v)-1)){ 
    ind1 <- (yhat>v[j]&yhat<=v[j+1])
    ind2 <- (s==0)
    ind3 <- (s==1)
    prob=(sum(ind1*ind3)/sum(ind3)*pp)/(sum(ind1*ind2)/sum(ind2))
    totaln=sum(ind1*ind2)
    value[j]=prob*totaln
  }
  return(value) 
}

para_ab<-function(s, predy, v, pp,cc1){
  yhat <- predy
  value <- numeric(length(v))
  for (j in 1:(length(v)-1)){ 
    ind1 <- (yhat>v[j]&yhat<=v[j+1])
    ind2 <- (s==0)
    num=sum(yhat*(1-cc1)*ind1)
    denom=sum((1-cc1*yhat)*ind1)
    totaln=sum(ind1*ind2)
    value[j]=num/denom*totaln
  }
  return(value) 
}

ppv<- function(y, predy, v){
  yhat <- predy
  ppv <- numeric(length(v))
  num_ppv <- numeric(length(v))
  denom_ppv <- numeric(length(v))
  for (j in 1:length(v)){
    ind <- (yhat>=v[j])
    num_ppv[j] <- sum(y[ind==1])
    denom_ppv[j] <-sum(ind)
    ppv[j]=num_ppv[j]/denom_ppv[j]
  }
  return(ppv)
}

npv<- function(y, predy, v){
  yhat <- predy
  npv <- numeric(length(v))
  num_npv <- numeric(length(v))
  denom_npv <- numeric(length(v))
  for (j in 1:length(v)){
    ind1 <- (yhat<v[j])
    ind2 <- (y==0)
    ind <- ind1*ind2
    num_npv[j] <- sum(ind)
    denom_npv[j] <-sum(ind1)
    npv[j]=num_npv[j]/denom_npv[j]
  }
  return(npv)
}

tpr<- function(y, predy, v){
  yhat <- predy
  tpr <- numeric(length(v))
  num_tpr <- numeric(length(v))
  denom_tpr <- numeric(length(v))
  for (j in 1:length(v)){
    ind1 <- (yhat>=v[j])
    ind2 <- (y==1)
    ind  <- ind1*ind2
    num_tpr[j] <- sum(ind)
    denom_tpr[j] <-sum(ind2)
    tpr[j]=num_tpr[j]/denom_tpr[j]
  }
  return(tpr)
}

fpr<- function(y, predy, v){
  yhat <- predy
  fpr <- numeric(length(v))
  num_fpr <- numeric(length(v))
  denom_fpr <- numeric(length(v))
  for (j in 1:length(v)){
    ind1 <- (yhat>=v[j])
    ind2 <- (y==0)
    ind  <- ind1*ind2
    num_fpr[j] <- sum(ind)
    denom_fpr[j] <-sum(ind2)
    fpr[j]=num_fpr[j]/denom_fpr[j]
  }
  return(fpr)
}
#######################################################
#Inference for accuracy measure function###############
#######################################################
asy_tpr_se <- function(ftrain, Betay,theta, cov, v, cutoff, delta ){

  tx=as.matrix(ftrain[,1:length(Betay)])
  nn=length(Betay)
  s=ftrain$s
  N=length(s)
  N_s1=sum(s)
  N_s0=N-N_s1
  
  beta=theta[1:nn]
  c=theta[nn+1]
  h=mean(s)
  pxb=plogis(tx %*% beta) #p(x;b)
  q=mean(pxb)
  p0=(q-h)/(1-h)
  
  ppxb=NULL #p'(x;b)
  for (i in 1:nn){
    ppxb=cbind(ppxb,pxb*(1-pxb)*tx[,i])
  }
  
  phi_b=NULL #n*nn
  for (i in 1:nn){
    phi_b=cbind(phi_b,s*ppxb[,i]/pxb-(1-s)*(c*ppxb[,i])/(1-c*pxb))
  }
  phi_cc=s*(1-c)-(1-s)*pxb/(1-c*pxb)*c*(1-c) 
  phi=t(cbind(phi_b,phi_cc)) 
  
  psi=N*cov%*%phi 
  psi_b=t(psi[1:nn,]) 
  psi_cc=t(psi[(nn+1),]) 

  IF=apply(psi,1,var)/N
  TC=diag(cov)
  V=diag(psi%*%t(psi)/N^2)
  
  #estimate accuracy measures in testing set
  ftest=ftrain
  ttx=as.matrix(ftest[,1:length(Betay)])
  tpxb=plogis(ttx %*% beta) 
  ts=ftest$s
  tN=length(ts)
  tN_s1=sum(ts)
  tN_s0=tN-tN_s1
  
  #estimation of accuracy measures
  #nonparametric version
  mm=length(v)
  tpr1=fpr1=ppv1=npv1=numeric(mm)
  tpr2=fpr2=ppv2=npv2=numeric(mm)
  psi_tpr1=psi_fpr1=psi_ppv1=psi_npv1=matrix(0,nrow=mm,ncol=length(s)) #v*nn
  psi_tpr2=psi_fpr2=psi_ppv2=psi_npv2=matrix(0,nrow=mm,ncol=length(s)) #v*nn
  
  yhat=tpxb
  ind2<- ts==0
  ind3<- ts==1
  for (j in 1:mm){ 
    ind1 <- (yhat>=v[j])
    ind4 <- 1-ind1
    pl0=mean(ind1*ind2)
    pl1=mean(ind1*ind3)
    npl0=mean(ind4*ind2)
    npl1=mean(ind4*ind3)
    
    tpr1[j]=pl1/mean(ind3)
    fpr1[j]=pl0/(1-q)-pl1*(q-h)/(h*(1-q))
    ppv1[j]=pl1/pl0*(q-h)/h
    npv1[j]=1-npl1/npl0*(q-h)/h
    
    psi_fpr1[j,]<-1/(1-q)*(ind1*ind2-pl0)-(q-h)/(h*(1-q))*(ind1*ind3-pl1)+pl1/(1-q)*q/h^2*(ind3-h)+
      ((pl0+pl1)/(1-q)^2-pl1/(h*(1-q)^2))*(yhat-q)
    psi_tpr1[j,]<-1/h*(ind1*ind3-pl1)-pl1/h^2*(ind3-h)
    psi_ppv1[j,]<-(q-h)/(pl0*h)*(ind1*ind3-pl1)-pl1*(q-h)/(pl0^2*h)*(ind1*ind2-pl0)+pl1/(pl0*h)*(yhat-q)-
      pl1/pl0*q/h^2*(ind3-h)
    psi_npv1[j,]<--(q-h)/(npl0*h)*(ind4*ind3-npl1)+npl1*(q-h)/(npl0^2*h)*(ind4*ind2-npl0)-npl1/(npl0*h)*(yhat-q)+
      npl1/npl0*q/h^2*(ind3-h)
    
  }
  
  auc1=NA
  auc2=NA
  
  #influence function of accuracy measures
  tpr_d1 =fpr_d1=ppv_d1=npv_d1=matrix(0,nrow=mm,ncol=nn) #v*nn
  tpr_d2 =fpr_d2=ppv_d2=npv_d2=matrix(0,nrow=mm,ncol=nn) #v*nn
  
  for (i in 1:nn){
    
    beta1=beta
    beta2=beta
    
    beta1[i]=beta1[i]+delta*beta[i]
    beta2[i]=beta2[i]-delta*beta[i]
    
    tpr_d1[,i] <- t((tpr_n(ts, ttx, beta1, v)-tpr_n(ts, ttx, beta2, v))/(2*delta*beta[i]))
    fpr_d1[,i] <- t((fpr_n(ts, ttx, beta1, v)-fpr_n(ts, ttx, beta2, v))/(2*delta*beta[i]))
    ppv_d1[,i] <- t((ppv_n(ts, ttx, beta1, v)-ppv_n(ts, ttx, beta2, v))/(2*delta*beta[i]))
    npv_d1[,i] <- t((npv_n(ts, ttx, beta1, v)-npv_n(ts, ttx, beta2, v))/(2*delta*beta[i])) 
  }
  
  tpr_nd= tpr_d1%*%t(psi_b)
  fpr_nd= fpr_d1%*%t(psi_b)
  ppv_nd= ppv_d1%*%t(psi_b)
  npv_nd= npv_d1%*%t(psi_b)
  
  tpr_an= psi_tpr1+tpr_nd
  fpr_an= psi_fpr1+fpr_nd
  ppv_an= psi_ppv1+ppv_nd
  npv_an= psi_npv1+npv_nd
  
  #asymptotic variance for non-para accuracy measures
  se_tpr1=sqrt(apply(tpr_an,1,var)/N)
  se_fpr1=sqrt(apply(fpr_an,1,var)/N)
  se_ppv1=sqrt(apply(ppv_an,1,var)/N)
  se_npv1=sqrt(apply(npv_an,1,var)/N)
  se_auc1=NA
  
  return(list( se_tpr1=se_tpr1,se_fpr1=se_fpr1,se_ppv1=se_ppv1,se_npv1= se_npv1,
               tpr1=tpr1,fpr1=fpr1,ppv1=ppv1,npv1=npv1))
  
}

psi_tpr_fpr<-function(beta,s,tx,psi_b, delta1, yhat, ts, q, h, v){
  mm=length(v)
  nn=length(beta)
  psi_tpr1=psi_fpr1=matrix(0,nrow=mm,ncol=length(ts)) #v*nn
  ind2<- ts==0
  ind3<- ts==1
  for (j in 1:mm){ 
    ind1 <- (yhat>=v[j])
    ind4 <- 1-ind1
    pl0=mean(ind1*ind2)
    pl1=mean(ind1*ind3)
    psi_fpr1[j,]<-1/(1-q)*(ind1*ind2-pl0)-(q-h)/(h*(1-q))*(ind1*ind3-pl1)+pl1/(1-q)*q/h^2*(ind3-h)+
      ((pl0+pl1)/(1-q)^2-pl1/(h*(1-q)^2))*(yhat-q)
    psi_tpr1[j,]<-1/h*(ind1*ind3-pl1)-pl1/h^2*(ind3-h)
  }
  
  tpr_d1 =fpr_d1=matrix(0,nrow=mm,ncol=nn) #v*nn
  
  for (i in 1:nn){
    
    beta1=beta
    beta2=beta
    
    beta1[i]=beta1[i]+delta1*beta[i]
    beta2[i]=beta2[i]-delta1*beta[i]
    
    #fpr, ppv, npv as function of beta and p0   
    tpr_d1[,i] <- t((tpr_n(s, tx, beta1, v)-tpr_n(s, tx, beta2, v))/(2*delta1*beta[i]))
    fpr_d1[,i] <- t((fpr_n(s, tx, beta1, v)-fpr_n(s, tx, beta2, v))/(2*delta1*beta[i]))
  }
  
  tpr_nd= tpr_d1%*%t(psi_b) #mm*n
  fpr_nd= fpr_d1%*%t(psi_b) #mm*n
  
  psi_tpr= psi_tpr1+tpr_nd
  psi_fpr= psi_fpr1+fpr_nd 
  
  return(list(psi_tpr=psi_tpr,psi_fpr=psi_fpr))
}  


psi_tpr_fpr_v<-function(beta,s,tx,psi_b, delta1, yhat, ts, q, h, v){
  nn=length(beta)
  ind2<- ts==0
  ind3<- ts==1
  ind1 <- (yhat>=v)
  ind4 <- 1-ind1
  pl0=mean(ind1*ind2)
  pl1=mean(ind1*ind3)
  psi_fpr1<-1/(1-q)*(ind1*ind2-pl0)-(q-h)/(h*(1-q))*(ind1*ind3-pl1)+pl1/(1-q)*q/h^2*(ind3-h)+
    ((pl0+pl1)/(1-q)^2-pl1/(h*(1-q)^2))*(yhat-q)
  psi_tpr1<-1/h*(ind1*ind3-pl1)-pl1/h^2*(ind3-h)
  
  tpr_d1 =fpr_d1=numeric(nn) 
  
  for (i in 1:nn){
    beta1=beta
    beta2=beta
    
    beta1[i]=beta1[i]+delta1*beta[i]
    beta2[i]=beta2[i]-delta1*beta[i]
    
    #fpr, ppv, npv as function of beta and p0   
    tpr_d1[i] <- t((tpr_n(s, tx, beta1, v)-tpr_n(s, tx, beta2, v))/(2*delta1*beta[i]))
    fpr_d1[i] <- t((fpr_n(s, tx, beta1, v)-fpr_n(s, tx, beta2, v))/(2*delta1*beta[i]))
  }
  
  tpr_nd= tpr_d1%*%t(psi_b) #n
  fpr_nd= fpr_d1%*%t(psi_b) #n
  
  psi_tpr= psi_tpr1+t(tpr_nd) #n*1
  psi_fpr= psi_fpr1+t(fpr_nd) #n*1
  
  return(list(psi_tpr=psi_tpr,psi_fpr=psi_fpr))
}

auc_se <- function(ftrain, Betay,theta, cov, cutoff, delta1, delta2){

  tx=as.matrix(ftrain[,1:length(Betay)])
  nn=length(Betay)
  s=ftrain$s
  N=length(s)
  N_s1=sum(s)
  N_s0=N-N_s1
  
  beta=theta[1:nn]
  c=theta[nn+1]
  h=mean(s)
  pxb=plogis(tx %*% beta) #p(x;b)
  q=mean(pxb)
  p0=(q-h)/(1-h)
  
  ppxb=NULL #p'(x;b)
  for (i in 1:nn){
    ppxb=cbind(ppxb,pxb*(1-pxb)*tx[,i])
  }
  
  phi_b=NULL 
  for (i in 1:nn){
    phi_b=cbind(phi_b,s*ppxb[,i]/pxb-(1-s)*(c*ppxb[,i])/(1-c*pxb))
  }
  phi_cc=s*(1-c)-(1-s)*pxb/(1-c*pxb)*c*(1-c) 
  phi=t(cbind(phi_b,phi_cc)) 
  
  psi=N*cov%*%phi 
  psi_b=t(psi[1:nn,]) 
  psi_cc=t(psi[(nn+1),]) 

  psi_tprfpr=psi_tpr_fpr(beta, s,tx,psi_b, delta1, pxb, s, q, h, cutoff)
  
  psi_tpr=psi_tprfpr$psi_tpr 
  psi_fpr=psi_tprfpr$psi_fpr 
  
  tpr1=tpr_n(s, tx, beta, cutoff)
  fpr1=fpr_n(s, tx, beta, cutoff)
  auc1=auc_self(fpr1,tpr1) 
  
  dfpr1=numeric(length(cutoff))
  dpsi_fpr1=matrix(0,nrow=length(cutoff),ncol=N)
  
  for (j in 1:length(cutoff)){
    v1=cutoff[j]+delta2
    v2=cutoff[j]
    #calculate dfpr/dv
    dfpr1[j]=(fpr_n(s, tx, beta, v1)-fpr_n(s, tx, beta, v2))/(delta2)
    
    #calculate d(psi_fpr)/dv  
    psi_fpr1=psi_tpr_fpr_v(beta,s,tx,psi_b, delta1, pxb, s, q, h, v1)$psi_fpr
    psi_fpr2=psi_tpr_fpr_v(beta,s,tx,psi_b, delta1, pxb, s, q, h, v2)$psi_fpr
    dpsi_fpr1[j,]=(psi_fpr1-psi_fpr2)/(delta2)
  }
  
  #calculate influence function for auc
  psi_auc1=numeric(nrow(tx))    
  
  for (i in 1:nrow(tx)) {
    psi_auc1[i]=mean(psi_tpr[,i]*dfpr1)+mean(dpsi_fpr1[,i]*tpr1)  
  }
  
  #asymptotic variance for non-para accuracy measures
  se_auc1=sqrt(var(psi_auc1)/N)
  
  return(list(se_auc1=se_auc1,auc1=auc1))
}

