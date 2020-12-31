library(Rsolnp)
library(nloptr)
library(nleqslv)
library(maxLik)
library(mvtnorm)
library(gmm)
library(pROC)
library("e1071")
library("MASS")
library(maxLik)
library(boot)

source("functions for calibration paper.R")
source("chi2 statistic_qhat.R")
source("Calibration statistic replication function_qhat.R")

n=300000 #population
nt=50000 #train
nv=500  #validation
ntest=50000 #test
nrep=1000 #replication
Beta4=c(5.4,0.2,-1.0,-2,0.4,-1.4,2.4,0.6,1.8,2.8) # Beta for prevalence of 20%
tol=1e-8
iter=1000
seed1=1234
v=c(0.2,0.3,0.4,0.5) #cutoff for accuracy measures
cutoff=seq(0,1,0.001) #cutoff for estimation of AUC
c=0.5 #p(s=1|y=1)
psd=0 #q is known

#cutoff for different K
vv1=seq(0,1,0.5) #k=2
vv2=seq(0,1,length.out = 6) #k=5
vv3=seq(0,1,length.out = 11) #k=10
vv4=seq(0,1,length.out = 13) #k=12
vv5=seq(0,1,length.out = 16) #k=15

#model misspecification
mc1=c(2,5,8)#miss weak predictors
mc2=c(3,6,9)#miss moderate predictors
mc3=c(2,3,5,6,8,9)#miss weak & moderate predictors

#estimate type I error rate or power of detecting misspecification

fit20=replik_calib_qhat(n, nt, nv,ntest, nrep, Betay=Beta4, random=T, c, cbeta, iter, tol, seed1, alpha=Beta4, vv=vv2, delta=0.1,mis=F,mc=mc1)
mean(fit20$mtstat)
mean(fit20$mpvalue<0.05)
fit20$qc
fit20$casetable

fit21=replik_calib_qhat(n, nt, nv,ntest, nrep, Betay=Beta4, random=T, c, cbeta, iter, tol, seed1, alpha=Beta4, vv=vv2, delta=0.1,mis=T,mc=mc1)
mean(fit21$mtstat)
mean(fit21$mpvalue<0.05)
fit21$qc
fit21$casetable

fit22=replik_calib_qhat(n, nt, nv,ntest, nrep, Betay=Beta4, random=T, c, cbeta, iter, tol, seed1, alpha=Beta4, vv=vv2, delta=0.1,mis=T,mc=mc2)
mean(fit22$mtstat)
mean(fit22$mpvalue<0.05)
fit22$qc
fit22$casetable

fit23=replik_calib_qhat(n, nt, nv,ntest, nrep, Betay=Beta4, random=T, c, cbeta, iter, tol, seed1, alpha=Beta4, vv=vv2, delta=0.1,mis=T,mc=mc3)

mean(fit23$mtstat)
mean(fit23$mpvalue<0.05)
fit23$casetable
fit23$qc


