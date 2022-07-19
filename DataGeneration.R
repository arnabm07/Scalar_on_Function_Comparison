library(mgcv); library(lme4); library(RLRsim); 
library(fda); library(fda.usc); library(sde)

########################################################################
########################################################################
##                                                                    ##
##  :::INPUTS:::                                                      ##
##                                                                    ##
##     n (or N): sample size                                          ##
##     m: number of grid points                                       ##
##     delta: parameter controlling the departure from the null       ##
##                                                                    ##
##  :::OUTPUTS:::                                                     ##
##                                                                    ##
##     Yi: scalar response                                            ##
##     Xit: functional covariate measured with or without error       ##
########################################################################
########################################################################

rm(list=ls())

genA<-function(n,m,delta,seed,norm=TRUE){   # norm=TRUE for including norm
  
  set.seed(112233+seed)
  t<-seq(0,1,len=m)
  
  #Generate functional covariate from "Ornstein Uhlenbeck process"
  Xit<-rproc2fdata(n, t=seq(0,1,len=m), mu=rep(0,length(t)), sigma="OU")
  tt<-fdata(mdata=rep(1,m),argvals=seq(0,1,len=m),rangeval=c(0,1))
  
  beta1t<-fdata(mdata=sin(2*pi*t)-cos(2*pi*t),argvals=seq(0,1,len=m),rangeval=c(0,1))
  if(norm) {Yi1<-inprod.fdata(Xit,beta1t)+delta*as.numeric(norm.fdata(Xit))+rnorm(n,mean=0,sd=0.1)}
  else {Yi1<-inprod.fdata(Xit,beta1t)+delta*inprod.fdata(Xit^2,tt)+rnorm(n,mean=0,sd=0.1)}
  
  Xit<-Xit$data 
  results<-list(Yi=Yi1, Xit=Xit) 
}


################################################################
################################################################
################################################################
################################################################
################################################################

rm(list=ls())

genB<-function(n,m,delta,seed,norm=TRUE){   # norm=TRUE for including norm
  
  set.seed(1122335+seed)
  t=seq(0,1,len=m)
  
  #Generate functional covariate from "Ornstein Uhlenbeck process"
  Xit<-rproc2fdata(n, t=seq(0,1,len=m), mu=rep(0,length(t)), sigma="OU")
  tt<-fdata(mdata=rep(1,m),argvals=seq(0,1,len=m),rangeval=c(0,1))
  
  beta2t<-fdata(mdata=t-(t-0.75)^2,argvals=seq(0,1,len=m),rangeval=c(0,1))
  if(norm) {Yi2<-inprod.fdata(Xit,beta2t)+delta*as.numeric(norm.fdata(Xit))+rnorm(n,mean=0,sd=0.1)}
  else {Yi2<-inprod.fdata(Xit,beta2t)+delta*inprod.fdata(Xit^2,tt)+rnorm(n,mean=0,sd=0.1)}
  
  Xit<-Xit$data   ##Extract Xit (a MATRIX)
  results<-list(Yi=Yi2, Xit=Xit)
}

################################################################
################################################################
################################################################
################################################################
################################################################

rm(list=ls())

genC<-function(N,delta,seed){
  
  set.seed(112233+seed)
  
  M=100		#	M is the number of points at which each curve is defined.
  
  Xmean=7		#	'Xmean' is the mean of X if simulated values are to be used.
  Ymean=4		#	'Ymean' is a paremeter but NOT necessarily the mean of Y.
  
  #	E will be the error terms.
  E = matrix(rnorm(N),N,1)
  
  X=matrix(rep(0,N*M),N,M)
  for (index in 1:N)
  {
    X[index,] = BM(x=Xmean, t0=0, T=1, M-1)
  }
  rm(index)
  
  #Xc<-X
  tij<-seq(0,1,length.out=M)
  #	'Xc' is the centered X.
  Xc = t(t(X) - colMeans(X))
  
  k=function(t)
  {
    output=t*0+1
    output
  }
  
  h=function(s,t,delta)
  {
    output = delta
    output
  }
  
  # h_matrix is a matrix of the values of h(.,.) at an M by M grid of points.
  
  h_matrix = matrix(0,M,M)
  for (row_index in 1:M)
  {
    for (col_index in 1:M)
    {
      h_matrix[row_index,col_index] = h(row_index/M,col_index/M,delta)
    }
  }
  rm(row_index,col_index)
  
  #	Y will be created from X and E
  Y=matrix(rep(0,N),N,1)
  for (index in 1:N)
  {
    Y[index] = Ymean + sum(k((1:M)/M)*Xc[index,])/M + 1/M^2 * t(Xc[index,]) %*% h_matrix %*% (Xc[index,]) + E[index] 
  }
  
  Y<-as.vector(Y)
  
  return(list(Yi=Y,Xit=Xc))
}

################################################################
################################################################
################################################################
################################################################
################################################################

##############################################################################################
##############################################################################################
#####  NOTE:"lambda" denotes "phi" in the paper.                                         #####     
#####                                                                                    #####     
#####  In our paper, we may need to denote "delta=1-phi" to keep the increasing order    #####
#####  when we present power results.                                                    #####
##############################################################################################
##############################################################################################
rm(list=ls())

source('CreateData.R')
source('TestFuns.R')

seed <- 112233
sig2e <- 1 # error variance
sigmax <- 0 # measurement error variance
twoVC <- TRUE
#snre <- 4
J <- 30 # number of observations used to generate data
Ji <- J # number of observations per curve
nxbf <- ntbf <- 10 # number of basis functions for x and t axes
ns <- 500 # sample size
Xtype <- "McLeanEtAl"
trueF <- c("linear2", "cos")  # true surface is convex combination of plane and cosine surface
method <- "Bonferroni"  # "simple", "twoVC", or "Bonferroni" 
oracle <- FALSE
removeCon <- FALSE
fun <- switch(method, simple=TestLinearitySimple, twoVC=TestKnowEqualVC, Bonferroni=testfun)
nc <- switch(method, simple=4, twoVC=4, 5)
bs <- "ps"
args <- list(family=gaussian(), 
             splinepars=list(k=c(nxbf,ntbf),m=list(c(2,2),c(2,2)),extendXrange=.01, bs=bs),
             REML=TRUE, tvals=seq(0, 1, l=J))
args$removeCon <- switch(method, simple=removeCon)  # removeCon not arg to
if (method != "simple" && bs == "ps")
  args$psplines <- TRUE
if (method == "twoVC")
  args$no.nuisvc.test.if.0est <- TRUE
cons <- c(1,10)
sx <- c(4,4)
st <- c(.5,.5)

nsim <- 5000
data<-list()
lambda <- 1 

for(k in 1:length(ns)){
  if(oracle){
    L <- matrix(1/J, ns[k], J)
    tmat <- matrix(args$tvals, ns[k], J, byrow=TRUE)
  }
  for(j in 1:length(lambda)){
    
    for(i in 1:nsim){
      
      datasim <- CreateDatat2cc(n=ns[k],nfine=J,tfine=NULL,trueF=trueF,seed=seed+nsim*3*(k-1)+nsim*(j-1)+i-1,
                                sigmae=sig2e,snre=NULL,sigmax=sigmax,snrx=NULL,Xtype=Xtype,
                                extendXrange=.01,propMissing=NULL,Ji=J,sx=sx,st=st,cons=cons,J=2,
                                ConvexCombo=TRUE,lambda=lambda[j])
      Yi <- datasim$y
      Xit <- datasim$Xorig
      data[[i]] <- list(Yi,Xit)
      names(data[[i]])<-c("Yi","Xit")
      print(i)
    }
  }
}

#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################

rm(list=ls())

genE<-function(n,delta,seed){
  
  set.seed(112233+seed)
  
  ksi1<-as.vector(rnorm(n,0,sqrt(4)))       #true scores_1
  ksi2<-as.vector(rnorm(n,0,sqrt(1)))       #true scores_2 
  
  tij<-seq(0,10,length.out=30)
  t<-do.call(rbind,replicate(n,t(tij),simplify=FALSE))
  mut<-t+sin(t)
  phi1t<-as.vector(-(1/sqrt(5))*cos(pi*tij/10))
  phi2t<-as.vector((1/sqrt(5))*sin(pi*tij/10))
  
  eps<-matrix(rnorm(n*length(tij),0,0.5), nrow=n)
  W<-as.matrix(mut + ksi1%*%t(phi1t) + ksi2%*%t(phi2t) + eps)
  
  Yi<-(ksi1+ksi2)+delta*(ksi1^2+ksi2^2+ksi1*ksi2)+sqrt(0.1)*rnorm(n)
  
  return(list(Yi=Yi,Xit=W))
}

#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################

genF<-function(n,m,delta,seed,norm=TRUE){   # norm=TRUE for including norm
  
  set.seed(112233+seed)
  t<-seq(0,1,len=m)
  
  #Generate functional covariate from "Ornstein Uhlenbeck process"
  Xit<-rproc2fdata(n, t=seq(0,1,len=m), mu=rep(0,length(t)), sigma="OU")
  tt<-fdata(mdata=rep(1,m),argvals=seq(0,1,len=m),rangeval=c(0,1))
  
  beta1t<-fdata(mdata=sin(2*pi*t)-cos(2*pi*t),argvals=seq(0,1,len=m),rangeval=c(0,1))
  Yi<-delta*inprod.fdata(Xit,beta1t)+rnorm(n,mean=0,sd=0.1)
  
  Xit<-Xit$data 
  results<-list(Yi=Yi, Xit=Xit) 
}

#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################

genG<-function(n,m,delta,seed){
  set.seed(112233+seed)
  
  Yi<-matrix(0,n)
  e1<-as.vector(rnorm(n,0,sd=8))           #true scores_1
  e2<-as.vector(rnorm(n,0,sd=2))           #true scores_2 
  e3<-as.vector(rnorm(n,0,sd=8/9))         #true scores_3
  e4<-as.vector(rnorm(n,0,sd=0.5))         #true scores_4
  
  t<-seq(0,1,len=m)
  
  phi1t<-as.vector(sin(pi*t))
  phi2t<-as.vector(cos(pi*t))
  phi3t<-as.vector(sin(2*pi*t))
  phi4t<-as.vector(cos(2*pi*t))
  
  Xit<-e1%*%t(phi1t)+e2%*%t(phi2t)+e3%*%t(phi3t)+e4%*%t(phi4t)
  
  F1<-function(x,t){2*x*sin(pi*t)}
  # F2<-function(x,t){10*cos(-0.125*x+0.25*t-5)}
  integ<-function(X,tt, FF){rowMeans(FF(X, tt))}
  
  tt<-matrix(t, nrow=n, ncol=m, byrow=TRUE)
  error<-rnorm(n,0,sqrt(1))
  Int_Y<-delta*integ(Xit, tt, F1)
  #Int_Y<-delta*integ(Xit, tt, F1) + (1-delta)*integ(Xit, tt, F2)
  Yi<-Int_Y + error   
  
  results<-list(Yi = Yi, Xit = Xit)  
  
}