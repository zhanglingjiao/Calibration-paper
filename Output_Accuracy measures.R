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
source("Accuracy replication function.R")

n=200000 #population
nt=80000 #train
nv=500 #validation
ntest=80000 #test
nrep=1000 #replication
Beta4=c(5.4,0.2,-1.0,-2,0.4,-1.4,2.4,0.6,1.8,2.8) #beta for 20% prevalence
c=0.5 #p(s=1|y=1)
cbeta=c(2,1.2,-.5,1,2,-.4,.4,1.6,1,.8)#not used
tol=1e-8
iter=1000
seed1=1234
v=seq(0.1,0.9,0.1) #cutoffs for accuracy measures
vv=seq(0,1,0.1) #cutoffs for computing para & non-para cases 
cutoff=seq(0,1,0.01) #used for auc estimation
mc1=c(2,5,8)#missing weak predictors, not used
mc2=c(2,3,5,6,8,9)#missing weak and moderate predictors, not used
c0=0.2 #only used for stratified scenarios, not used here
c1=0.8 #only used for stratified scenarios, not used here

##with bootstrap, get estimate of beta
fit=replik_accu(n, nt, nv,ntest, nrep, Beta4, random=TRUE, c, cbeta, iter, tol, seed1, method='LZ_c_optim', v, mis=FALSE, mc,
               stratify=FALSE, c0, c1, cutoff,vv,delta1=0.1,delta2=0.01)

fit$estimation #parameter estimation

fit$accu_strain_ase1 #asymptotic standard error for self-defined accuracy measures
fit$accu_ettrain_sd #empirical standard error for benchmark accuracy measures

fit$accu_strain #self-defined accuracy measures
fit$accu_ettrain #benchmark accuracy measures
fit$accu_auc #AUC (SE)



