
replik_calib_slope=function(n, nt, nv,ntest, nrep, Betay, random, c, cbeta, iter, tol, seed1, alpha, vv, delta,mis,mc){
  
  set.seed(seed1)
  
  slope1=slope2=list()
  
  tstat= total=pvalue=list()
  tU=tV=list()
  tUn=tUm=list()
  pcase=npcase=list()
  qq=cc=list()

  replication=0
  
  for (i in 1:nrep){
    
    name <- paste('rep',i,sep='')
    replication=replication+1

    if(i==200){print('replication=200')}
    if(i==400){print('replication=400')}
    if(i==600){print('replication=600')}
    if(i==800){print('replication=800')}
    conv1=0
    
    while (conv1==0){   
    sample=rsample_ns(n,nt,nv,ntest,Betay,random,c,cbeta)
    prevalence=sample$prevalence #true prevalence
    
    test=sample$test
    valid=sample$valid
    train=sample$train
  
    if (mis==FALSE) {        
        mBetay=Betay
    }
      
    if (mis==TRUE) { 
          mBetay=Betay[-mc]
          train=train[,-mc]
          test=test[,-mc]
    }    

       train1=as.matrix(train[,1:length(mBetay)])
       train$lproby=train1%*% mBetay

       #construct positive only data, (inter,logit_p,s)       
       train_po=data.frame(cbind(rep(1,nrow(train)),train$lproby,train$s))
       names(train_po)=c('inter','lphat','s')
  
       nBetay=c(0.5,0.5)
       fit_po=tryCatch(LZ_c_optim(train_po, iter, tol, nBetay), error=function(err) list(conv=0))
       conv1=fit_po$conv
    }

       slope1[[name]]=fit_po$theta[1:length(nBetay)]
       
       ##slope, logit y~proby
       ylogit <- glm(y ~ lproby, data = train, family = "binomial")
       slope2[[name]]= coef(ylogit)  
       
}
  #estimation 
  mslope1=do.call("rbind",slope1)
  
  slope1_m=apply(mslope1,2,mean)
  slope1_sd=apply(mslope1,2,sd)
  
  mslope2=do.call("rbind",slope2)
  
  slope2_m=apply(mslope2,2,mean)
  slope2_sd=apply(mslope2,2,sd)

  summary=cbind(slope1_m,slope2_m,slope1_sd,slope2_sd)
  colnames(summary)=c('ml_m','y_m','ml_sd','y_sd')
  
  #output
  return(list(slope1_m=slope1_m,slope2_m=slope2_m,summary=summary))
}

