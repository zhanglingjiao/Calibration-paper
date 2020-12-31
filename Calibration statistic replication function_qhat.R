
replik_calib_qhat=function(n, nt, nv,ntest, nrep, Betay, random, c, cbeta, iter, tol, seed1, alpha, vv, delta,mis,mc){
  
  set.seed(seed1)
  
  tstat= total=pvalue=list()
  tU=tV=list()
  tUn=tUm=list()
  pcase=npcase=list()
  qq=cc=list()

  replication=0
  
  for (i in 1:nrep){
    
    name <- paste('rep',i,sep='')
    replication=replication+1
    conv1=0
    
    if(i==200){print('replication=200')}
    if(i==400){print('replication=400')}
    if(i==600){print('replication=600')}
    if(i==800){print('replication=800')}
    
    sample=rsample_ns(n,nt,nv,ntest,Betay,random,c,cbeta)
    prevalence=sample$prevalence #true prevalence
    
    test=sample$test
    valid=sample$valid
    train=sample$train

  while (conv1==0){
        #sample presence, background and test data  
    if (mis==FALSE) {        
        fit=tryCatch(LZ_c_optim(train, iter, tol,Betay), error=function(err) list(conv=0))
        conv1=fit$conv
        mBetay=Betay
    }
      
    if (mis==TRUE) { 
          train=train[,-mc]
          test=test[,-mc]
          mBetay=Betay[-mc]
          #Fit model
          fit=tryCatch(LZ_c_optim(train, iter, tol,mBetay), error=function(err) list(conv=0))
          conv1=fit$conv
    }    
}
       #extract estimates from optim
       coef1=fit$theta
       cov1=fit$cov
       bcov1=fit$bcov
       bbeta=as.matrix(coef1[1:length(mBetay)])
       
       ##prepare accuracy measures in training set
       train1=as.matrix(train[,1:length(mBetay)])
       train$proby=plogis(train1%*% bbeta)

    #compute estimated prevalence using training set only
        q1=mean(train$proby)
        cc1=coef1[length(mBetay)+1]

    #calibration
    calibration0=calib_qhat(train, mBetay, coef1, cov1, bbeta, vv, delta)
    stat0=calibration0$stat
    V0=calibration0$V
    U0=calibration0$U
    Un0=calibration0$Un
    Um0=calibration0$Um
    vtotal0=calibration0$vtotal
    pvalue0=calibration0$pvalue
    pcase0=calibration0$case_para
    npcase0=calibration0$case_npara
    
    tstat[[name]] <- stat0
    tV[[name]] <- V0
    tU[[name]] <- U0
    tUn[[name]] <- Un0
    tUm[[name]] <- Um0
    total[[name]] <- vtotal0
    pvalue[[name]] <- pvalue0
    pcase[[name]] <- pcase0
    npcase[[name]] <- npcase0
    qq[[name]]<-q1
    cc[[name]]<-cc1
  
}
  #estimation 
  mtstat=do.call("rbind",tstat)
  mtU=do.call("rbind",tU)
  mtUn=do.call("rbind",tUn)
  mtUm=do.call("rbind",tUm)
  mtotal=do.call("rbind",total)
  mpvalue=do.call("rbind",pvalue)
  mpcase=do.call("rbind",pcase)
  mnpcase=do.call("rbind",npcase)
  mqq=do.call("rbind",qq)
  mcc=do.call("rbind",cc)
  
  
  qq_m=mean(mqq)
  qq_se=sd(mqq)
  
  cc_m=mean(mcc)
  cc_se=sd(mcc)
  
  qc=c(qq_m,qq_se,cc_m,cc_se)
  names(qc)=c('q_m','q_se','c_m','c_se')
  
  case_p=apply(mpcase,2,mean)
  case_np=apply(mnpcase,2,mean)
  
  casetable=cbind(case_np,case_p) 
  
 
if (length(vv)>3){
  mtv=apply(simplify2array(tV), 1:2, mean)
}
  if (length(vv)==3){
  mtv=mean(do.call("rbind",tV))
}  
  #output
  return(list(mtstat=mtstat,mtU=mtU, tV=tV, mtv=mtv,casetable=casetable,qc=qc,
              mtUn=mtUn, mtUm=mtUm, mtotal=mtotal, mpvalue=mpvalue, mpcase=mpcase, mnpcase=mnpcase))
}

