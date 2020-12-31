
calib_qhat=function(ftrain, Betay, theta, cov, alpha, v, delta ){
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
  
  ppxb=NULL 
  for (i in 1:nn){
    ppxb=cbind(ppxb,pxb*(1-pxb)*tx[,i])
  }
  
  phi_b=NULL 
  for (i in 1:nn){
    phi_b=cbind(phi_b,s*ppxb[,i]/pxb-(1-s)*(c*ppxb[,i])/(1-c*pxb))
  }
  phi_cc=s*(1-c)-(1-s)*pxb/(1-c*pxb)*c*(1-c) #n*1
  
  phi=t(cbind(phi_b,phi_cc)) 
  
  psi=N*cov%*%phi 
  psi_b=t(psi[1:nn,]) 

  gx<-plogis(tx %*% alpha)

  mm=length(v)-1
  npara=mpara=numeric(mm)
  phi_an=phi_am=matrix(0,nrow=mm,ncol=N)
 
  ind2 <- (s==0)
  ind3 <- (s==1)
  vtotal=NULL
  vpl0=vpl1=vpl=NULL

for (j in 1:mm){ 
    ind1 <- (gx>v[j]&gx<=v[j+1])
    total=sum(ind1*ind2)
    vtotal=c(vtotal,total)
    
    npara[j]=mean(ind1*ind3)*(q-h)/(h*mean(ind1*ind2))
    mpara[j]=sum(pxb*(1-c)*ind1)/sum((1-c*pxb)*ind1)

    pl0=mean(ind1*ind2)
    pl1=mean(ind1*ind3)
    pl=mean(ind1)

    psi_n0=(q-h)/(h*pl0)*(ind1*ind3-pl1)-
          (q-h)*pl1/(h*pl0^2)*(ind1*ind2-pl0)-
          (q*pl1)/(h^2*pl0)*(ind3-h)+
          pl1/(h*pl0)*(pxb-q)
    psi_m0=((pxb*ind1*(1-c))-mpara[j]*(1-c*pxb)*ind1)/mean((1-c*pxb)*ind1)

    npara_d=mpara_d=numeric(nn)

    for (i in 1:nn){
      beta1=beta
      beta2=beta
      
      beta1[i]=beta[i]+delta
      beta2[i]=beta[i]-delta
      
      tpxb1=plogis(tx%*% beta1)
      tpxb2=plogis(tx%*% beta2)
      
      q_1=mean(tpxb1)
      q_2=mean(tpxb2)
     
      npara_1=(mean(ind1*ind3)*(q_1-h))/(h*mean(ind1*ind2))
      npara_2=(mean(ind1*ind3)*(q_2-h))/(h*mean(ind1*ind2))
     
      mpara_1=sum(tpxb1*(1-c)*ind1)/sum((1-c*tpxb1)*ind1) 
      mpara_2=sum(tpxb2*(1-c)*ind1)/sum((1-c*tpxb2)*ind1) 
     
      npara_d[i] <- (npara_1-npara_2)/(2*delta)
      mpara_d[i] <- (mpara_1-mpara_2)/(2*delta)
      
    }   
     c1=c+delta
     c2=c-delta
     mpara_c1=sum(pxb*(1-c1)*ind1)/sum((1-c1*pxb)*ind1)
     mpara_c2=sum(pxb*(1-c2)*ind1)/sum((1-c2*pxb)*ind1)
     mpara_c <- (mpara_c1-mpara_c2)/(2*delta)
     mpara_cc=mpara_c*c*(1-c)
    
    mpara_theta=c(mpara_d,mpara_cc)
    
    a_m=t(psi)%*%mpara_theta
    a_n= psi_b%*%npara_d

    phi_an[j,]=psi_n0+a_n 
    phi_am[j,]=psi_m0+a_m
}

Un=npara[-mm]
Um=mpara[-mm]
U=Un-Um

case_para=mpara*vtotal
case_npara=npara*vtotal

nphi_an=phi_an[-mm,]
nphi_am=phi_am[-mm,]
nphi_nm=nphi_an-nphi_am

if (mm==2){
V=var(nphi_nm)/N
stat=U^2/V
}

if (mm>2){
V=nphi_nm%*%t(nphi_nm)/N^2
stat=t(U)%*%solve(V)%*%U
}

pvalue=pchisq(stat, df=mm-1, lower.tail=FALSE)
return(list(stat=stat,U=U,Un=Un,Um=Um,V=V,vtotal=vtotal,pvalue=pvalue,case_para=case_para,case_npara=case_npara))
}

