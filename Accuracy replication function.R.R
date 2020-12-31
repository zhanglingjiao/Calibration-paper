
replik_accu=function(n, nt, nv,ntest, nrep, Betay, random, c, cbeta, iter, tol, seed1, method, v, mis, mc,
                stratify, c1, c2, cutoff,vv,delta1,delta2){
  
  set.seed(seed1)
  
  scoef=sse=sconv=list()
  sc=asec=list()
  sq=aseq=list()

  tppv=ttpr=tnpv=tfpr=tauc=list()#non-parametric self-defined accuracy on training set
  etppv=ettpr=etnpv=etfpr=etauc=list()#empirical accuracy on training set
  tppv_ase1=ttpr_ase1=tnpv_ase1=tfpr_ase1=tauc_ase1=list()

  replication=0
  
  for (i in 1:nrep){
    name <- paste('rep',i,sep='')
    replication=replication+1
    conv1=0
    
    if(i==200){print('replication=200')}
    if(i==400){print('replication=400')}
    if(i==600){print('replication=600')}
    if(i==800){print('replication=800')}

while (conv1==0){
      sample=rsample_ns(n,nt,nv,ntest,Betay,random,c,cbeta)
      prevalence=sample$prevalence #true prevalence
      test=sample$test
      valid=sample$valid
      train=sample$train
  
      #Fit model
      if (method=='LZ_c_optim'){  
        fit=tryCatch(LZ_c_optim(train, iter, tol,Betay), error=function(err) list(conv=0))
      }
      if (method=='LZ_c2_optim'){  
        fit=tryCatch(LZ_c2_optim(train, iter, tol,Betay), error=function(err) list(conv=0))
      }
      conv1=fit$conv
    }
    #extract estimates from optim
    coef1=fit$theta
    aseq1=fit$qse
    se1=fit$se
    cov1=fit$cov
    bcov1=fit$bcov
    bbeta=as.matrix(coef1[1:length(Betay)])
    
    train1=as.matrix(train[,1:length(Betay)])
    train$proby=plogis(train1%*% bbeta)
    
    test1=as.matrix(test[,1:length(Betay)])  
    test$proby=plogis(test1%*% bbeta)
    
    ##restrict to s=0
    nntrain=train[train$s==0,]
    nntest=test[test$s==0,]
    
    q1=mean(train$proby)
    cc1=coef1[length(Betay)+1]
    sec1=se1[length(Betay)+1]
    seq1=aseq1
    pp=(q1-mean(train$s))/(1-mean(train$s))
    
    ##compute para & non-para case in S=0
    paracase0=para_ab(test$s, test$proby, vv, pp,cc1)
    nonparacase0=nonpara_ab(test$s, test$proby, vv, pp,cc1)

#se of tpr, fpr, ppv, npv
    fasy_ase1=asy_tpr_se(train, Betay,coef1, cov1, v, cutoff,delta1)
    
#se of auc
    fauc_ase1=auc_se(train, Betay,coef1, cov1, cutoff, delta1, delta2)

#Self-defined non-para accuracy measure
    tpr0=fasy_ase1$tpr1
    ppv0=fasy_ase1$ppv1
    fpr0=fasy_ase1$fpr1
    npv0=fasy_ase1$npv1
    auc0=fauc_ase1$auc1
#Empirical accuracy measure on s=0
    tpr03=tpr(nntrain$y, nntrain$proby, v)
    ppv03=ppv(nntrain$y, nntrain$proby, v)
    fpr03=fpr(nntrain$y, nntrain$proby, v)
    npv03=npv(nntrain$y, nntrain$proby, v)    
    auc03=auc_self(fpr(nntrain$y, nntrain$proby, cutoff),tpr(nntrain$y, nntrain$proby, cutoff)) 
#Asy se
    tpr_ase01=fasy_ase1$se_tpr1
    fpr_ase01=fasy_ase1$se_fpr1
    ppv_ase01=fasy_ase1$se_ppv1
    npv_ase01=fasy_ase1$se_npv1
    auc_ase01=fauc_ase1$se_auc1

    scoef[[name]] <- coef1
    sse[[name]] <- se1
    sq[[name]] <- q1
    sc[[name]] <- cc1
    aseq[[name]] <- seq1
    asec[[name]] <- sec1

    tppv[[name]] <- ppv0
    ttpr[[name]] <- tpr0
    tnpv[[name]] <- npv0
    tfpr[[name]] <- fpr0
    tauc[[name]] <- auc0

    tppv_ase1[[name]] <- ppv_ase01
    ttpr_ase1[[name]] <- tpr_ase01
    tnpv_ase1[[name]] <- npv_ase01
    tfpr_ase1[[name]] <- fpr_ase01
    tauc_ase1[[name]] <- auc_ase01

    etppv[[name]] <- ppv03
    ettpr[[name]] <- tpr03
    etnpv[[name]] <- npv03
    etfpr[[name]] <- fpr03
    etauc[[name]] <- auc03

  }
  scoef=do.call("rbind",scoef)
  sse=do.call("rbind",sse)
  sq=do.call("rbind",sq)
  sc=do.call("rbind",sc)
  aseq=do.call("rbind",aseq)
  asec=do.call("rbind",asec)
  
  coef=apply(na.omit(scoef),2,mean)
  ccoef=c(Betay,c) 
  bias=coef-ccoef

  ese=apply(na.omit(scoef),2,sd) #empirical se for coef
  mse=apply(na.omit(sse),2,mean) #asym se for coef
  
  estimation=data.frame(coef,bias,mse,ese)
  colnames(estimation) <- c('coef','bias','ase','ese')
  
  #mean estimate for q&c
  meanq=apply(na.omit(sq),2,mean) 
  meanc=apply(na.omit(sc),2,mean)  
  #empirical se for q &c
  eseq=apply(na.omit(sq),2,sd)
  esec=apply(na.omit(sc),2,sd)
  #asymptotic se for q &c
  maseq=apply(na.omit(aseq),2,mean)
  masec=apply(na.omit(asec),2,mean)
  
  cq_true=c(c,prevalence)
  cq_hat=c(meanc,meanq)
  cq_bias=cq_hat-cq_true
  cq_ase=c(masec,maseq)
  cq_ese=c(esec,eseq)
  
  cq_est=data.frame(cq_hat,cq_bias,cq_ase,cq_ese)
  colnames(cq_est) <- c('coef','bias','ase','ese')
  
#self-define accuracy measures  
  mtppv=do.call("rbind",tppv)
  mttpr=do.call("rbind",ttpr)
  mtnpv=do.call("rbind",tnpv)
  mtfpr=do.call("rbind",tfpr)
  mtauc=do.call("rbind",tauc)
  
  tppv_ase1=do.call("rbind",tppv_ase1)
  ttpr_ase1=do.call("rbind",ttpr_ase1)
  tnpv_ase1=do.call("rbind",tnpv_ase1)
  tfpr_ase1=do.call("rbind",tfpr_ase1)
  tauc_ase1=do.call("rbind",tauc_ase1)

#empirical accuracy measures  
  metppv=do.call("rbind",etppv)
  mettpr=do.call("rbind",ettpr)
  metnpv=do.call("rbind",etnpv)
  metfpr=do.call("rbind",etfpr)
  metauc=do.call("rbind",etauc)
  
  if (length(v)>1){

#empirical standard errors on self-defined accuracy measures
    mtppv_sd=apply(na.omit(mtppv),2,sd)
    mtnpv_sd=apply(na.omit(mtnpv),2,sd)
    mttpr_sd=apply(na.omit(mttpr),2,sd)
    mtfpr_sd=apply(na.omit(mtfpr),2,sd)
    mtauc_sd=sd(na.omit(mtauc))

#asymptotic standard errors on self-defined accuracy measures 
    mtppv_ase1=apply(na.omit(tppv_ase1),2,mean)
    mttpr_ase1=apply(na.omit(ttpr_ase1),2,mean)
    mtnpv_ase1=apply(na.omit(tnpv_ase1),2,mean)
    mtfpr_ase1=apply(na.omit(tfpr_ase1),2,mean) 
    mtauc_ase1=mean(na.omit(tauc_ase1))

#empirical standard error on true accuracy measures
    metppv_sd=apply(na.omit(metppv),2,sd)
    metnpv_sd=apply(na.omit(metnpv),2,sd)
    mettpr_sd=apply(na.omit(mettpr),2,sd)
    metfpr_sd=apply(na.omit(metfpr),2,sd)
    metauc_sd=sd(metauc)

#self-defined accuracy measures
    mtppv=apply(na.omit(mtppv),2,mean)
    mttpr=apply(na.omit(mttpr),2,mean)
    mtnpv=apply(na.omit(mtnpv),2,mean)
    mtfpr=apply(na.omit(mtfpr),2,mean)
    mtauc=mean(na.omit(mtauc))
    
#true accuracy measures   
    metppv=apply(na.omit(metppv),2,mean)
    mettpr=apply(na.omit(mettpr),2,mean)
    metnpv=apply(na.omit(metnpv),2,mean)
    metfpr=apply(na.omit(metfpr),2,mean)
    metauc=mean(metauc)
    }
####accuracy measures, estimation 
  #self-defined
  accu_strain=data.frame(v,mtppv,mtnpv,mttpr,mtfpr )
  colnames(accu_strain) <- c('cutoff','ppv','npv','tpr','fpr')

  #empirical  
  accu_ettrain=data.frame(v,metppv,metnpv,mettpr,metfpr)
  colnames(accu_ettrain) <- c('cutoff','ppv','npv','tpr','fpr')
 
#####standard errors 
  #self-defined, asymptotic standard error 
  accu_strain_ase1=data.frame(v,mtppv_ase1,mtnpv_ase1,mttpr_ase1,mtfpr_ase1)
  colnames(accu_strain_ase1) <- c('cutoff','ppv','npv','tpr','fpr')
  
  #self-defined, empirical standard error
  accu_strain_sd=data.frame(v,mtppv_sd,mtnpv_sd,mttpr_sd,mtfpr_sd)
  colnames(accu_strain_sd) <- c('cutoff','ppv','npv','tpr','fpr') 
  
  #benchmark (true), empirical
  accu_ettrain_sd=data.frame(v,metppv_sd,metnpv_sd,mettpr_sd,metfpr_sd)
  colnames(accu_ettrain_sd) <- c('cutoff','ppv','npv','tpr','fpr')
  
  #auc
  accu_auc=data.frame(mtauc,metauc,mtauc_ase1,mtauc_sd,metauc_sd)
  names(accu_auc) <- c('Est','Empirical','ASE','ESE','ESE_Emp')
  
  #output
  return(list(estimation=estimation,cq_est=cq_est,
              accu_strain=accu_strain, accu_ettrain=accu_ettrain, 
              accu_ettrain_sd=accu_ettrain_sd,accu_strain_sd=accu_strain_sd, accu_strain_ase1=accu_strain_ase1,accu_auc=accu_auc))
}

