library(fda.usc)

########################################################################
###
###  Script created by Merve T, Aug 1, 2016.
###
###  Used in the simulations, null test comp by MYT, MAC, AMS, AM.
########################################################################

GarciaEtal_null<-function(Yi,Xit){
  
  Xit<-fdata(mdata=Xit,argvals=seq(0,1,len=ncol(Xit)),rangeval=c(0,1))
  
  beta0<-rep(0,length(Xit$argvals))
  
  pcvm.sim<-flm.test(Xit,Yi,B=5000,
                     beta0.fdata=fdata(mdata=beta0,argvals=Xit$argvals,rangeval=Xit$rangeval),
                     type.basis = "bspline",
                     est.method="basis", plot.it=FALSE)
  p.val<-pcvm.sim$p.value
  results<-p.val
  return(results)
}


########################################################################
###
###  Script created by Janet Kim and Merve T, Aug 1, 2016.
###
###  Used in the simulations, null test comp by MYT, MAC, AMS, AM.
########################################################################

McLeanEtal_null<-function(Yi,Xit){
  
  n = length(Yi)
  m = ncol(Xit)
  Kx = Kt = 10
  t = seq(0,1, len=m)
  t.vec = rep(t, n)
  Xit.vec = as.vector(t(Xit))
  
  # Define boxproducts:
  boxprod = function(a,b){
    if (is.matrix(a) & is.matrix(b)){result = (a%x%matrix(1, nrow=1, ncol=ncol(b)))*(matrix(1, nrow=1, ncol=ncol(a))%x%b)}
    if (is.matrix(a) & is.vector(b)){result = (a%x%1)*(matrix(1, nrow=1, ncol=ncol(a))%x%b)}
    if (is.vector(a) & is.matrix(b)){result = (a%x%matrix(1, nrow=1, ncol=ncol(b)))*(1%x%b)}
    if (is.vector(a) & is.vector(b)){result = (a%x%1)*(1%x%b)}
    result
  }
  
  # Define weight matrix based on the Midpoint rule:
  i.vec<-rep(1,n)
  Iden<-diag(i.vec)
  ones<-rep(1,m)
  L.mat = (1/m)*(Iden%x%t(ones))
  
  # Construct tensor product smooth:
  Bx.obj = create.bspline.basis(rangeval=range(Xit), nbasis=Kx, norder=4)
  Bt.obj = create.bspline.basis(rangeval=range(t), nbasis=Kt, norder=4)
  Bx = eval.basis(Xit.vec, Bx.obj)             # n*m-by-Kx
  Bt = eval.basis(t.vec, Bt.obj)                # n*m-by-Kt
  BxBt = tensor.prod.model.matrix(list(Bx, Bt)) # n*m-by-Kx*Kt
  Px = eval.penalty(basisobj=Bx.obj, Lfdobj=int2Lfd(2))  #
  Pt = eval.penalty(basisobj=Bt.obj, Lfdobj=int2Lfd(2))  #
  
  # Eigendecomposition of marginal penalty:
  
  # Direction x:
  eig.x = eigen(Px)
  Upx = eig.x$vectors[,1:(Kx-2)]
  Unx = eig.x$vectors[,c(Kx-1, Kx)]
  D0.5inv.x = diag(1/sqrt(eig.x$values[1:(Kx-2)]))
  # Direction t:
  eig.t = eigen(Pt)
  Upt = eig.t$vectors[,1:(Kt-2)]
  Unt = eig.t$vectors[,c(Kt-1, Kt)]
  D0.5inv.t = diag(1/sqrt(eig.t$values[1:(Kt-2)]))
  
  Xx = cbind(1,Xit.vec)        # n*m-by-2             
  Zx = Bx%*%Upx%*%D0.5inv.x     # n*m-by-(Kx-2)
  Xt = cbind(1, t.vec)         # n*m-by-2
  Zt = Bt%*%Upt%*%D0.5inv.t     # n*m-by-(Kt-2)
  
  model.fix = boxprod(Xx,Xt)[,-2]   # [1:t:x:t*x]; n*m-by-3
  model.r1 = boxprod(Xit.vec, Zt)  # [x*Zt1:...:x*Zt8]; n*m-by-8
  model.r2 =  boxprod(Zx, Xt)       # n*m-by-8*2
  model.r3 = boxprod(Zx, Zt)        # n*m-by-8*8
  
  # multiply by wieght (L.mat)
  LXx = L.mat%*%model.fix       # n-by-3
  LZ1 = L.mat%*%model.r1        # n-by-8
  LZ2 = L.mat%*%model.r2        # n-by-8*2
  LZ3 = L.mat%*%model.r3        # n-by-8*8
  LZ  = cbind(LZ1, LZ2, LZ3)     # n-by-88
  
  dat<-data.frame(Yi,LXx,LZ1)
  index<-rep(1,n)
  U1<-LXx[,2]
  U2<-LXx[,3]
  colnames(dat)<-c("Yi","index","U1","U2",paste("z",1:(ncol(LZ1)),sep=""))
  random.model<-as.formula(paste("~-1+",paste(paste("z",1:(ncol(LZ1)),sep=""),collapse="+")))
  NULL1<-lm(Yi~1, data=dat)
  ALT1<-lme(Yi~1+U1+U2, random=list(index=pdIdent(random.model)), method="ML", data=dat)
  
  pval<-exactLRT(m=ALT1,m0=NULL1)$p
  results<-pval
  
  return(results)
  
}


########################################################################
###
###  Script created by Merve T, Aug 1, 2016.
###
###  Used in the simulations, null test comp by MYT, MAC, AMS, AM.
########################################################################

KSMEtal<-function(Yi,Xit){
  
  n<-length(Yi)
  m<-ncol(Xit)
  In<-diag(n)   #nxn identity matrix
  
  fpca<-fpca.sc(Xit,pve=0.99, var=TRUE)
  sn<-fpca$npc        #number of PCs
  M<-as.matrix(fpca$scores/sqrt(m))   #M matrix
  
  ones<-as.vector(rep(1,n))   #a vector with ones
  
  B<-cbind(ones,M)   #B matrix
  
  #Projection matrices
  PB<-B%*%chol2inv(chol(t(B)%*%B))%*%t(B)
  P1<-ones%*%t(ones)/n
  
  #sig.ML<-(t(Yi)%*%(In-P1)%*%Yi)/n   #MLE of sigma.square
  
  #Residual sum of squares under the alternative and null hyp.
  RSS.full<-t(Yi)%*%(In-PB)%*%Yi
  RSS.red<-t(Yi)%*%(In-P1)%*%Yi
  
  TF<-((RSS.red-RSS.full)/sn)/(RSS.full/(n-sn-1))   #F test statistic
  
  #####  OR  ####
  #TF<-(t(Yi)%*%(PB-P1)%*%Yi/sn)/(t(Yi)%*%(In-PB)%*%Yi/(n-sn-1)) 
  
  p.val<-1-pf(TF, sn, n-sn-1)
  
  results<-p.val
  return(results)
}
