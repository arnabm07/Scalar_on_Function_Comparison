########################################################################
###
###  Script created by Merve Tekbudak, Aug 1, 2016.
###
###  Used in the simulations, linear test comp by MYT, MAC, AMS, AM.
########################################################################

library(fda.usc)


GarciaEtal<-function(Yi,Xit){
  
  Xit<-fdata(mdata=Xit, argvals=seq(0,1,len=ncol(Xit)), rangeval=c(0,1))
  
  pcvm.sim<-flm.test(Xit, Yi, B=5000, beta0.fdata = NULL, type.basis = "bspline",
                     est.method="basis", plot.it=FALSE)
  p.val<-pcvm.sim$p.value
  results<-p.val
  return(results)
}


########################################################################
###
###  Script created by Janet Kim and Merve T, Aug 1, 2016.
###
###  Used in the simulations, linear test comp by MYT, MAC, AMS, AM.
########################################################################

McLeanEtal<-function(Yi,Xit){
  
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
  
  Z=list(LZ1, LZ2, LZ3)
  Sigma <- lapply(Z, function(z) diag(ncol(z)))
  
  require(RLRsim)
  
  block.diag <-
    function(...)
    {
      matrices <- list(...)
      k <- length(matrices)
      nrows <- sapply(matrices, nrow)
      ncols <- sapply(matrices, ncol)
      
      index.row <- rep(seq_len(k), nrows)
      index.col <- rep(seq_len(k), ncols)
      M <- matrix(0, sum(nrows), sum(ncols))
      for (i in seq_len(k)) {
        M[index.row == i, index.col == i] <- matrices[[i]]
      }
      
      M
    }
  
  m0 <- 1
  j0 <- seq_len(m0)
  j1 <- (m0 + 1L) : length(Z)
  m1 <- m0 + 1L
  
  Z[[m1]] <- do.call(cbind, Z[j1])
  Sigma[[m1]] <- do.call(block.diag, Sigma[j1])
  
  index<-rep(1,n)
  
  fit<-lme(Yi~-1+X, random=list(index=pdIdent(~-1+Z1), index=pdIdent(~-1+Z2)), method="REML")
  
  Z1 <- Z[[m1]]
  fit1 <- lme(Yi ~ -1 + LXx, random=list(index=pdIdent(~-1+LZ1), index=pdIdent(~-1+LZ2), index=pdIdent(~-1+LZ3)), method="REML")
  fit0 <- lme(Yi ~ -1 + LXx, random=list(index=pdIdent(~-1+LZ1)), method="REML")
  
  stat.obs <- 2 * (fit1$logLik - fit0$logLik)
  decomp <- eigen(Sigma[[m1]], TRUE)
  v <- decomp$values
  v[v < 0] <- 0
  
  nsim <- 10000
  seed <- 3142018
  
  R <- decomp$vectors %*% (sqrt(v) * t(decomp$vectors))
  stat.sim <- RLRsim::RLRTSim(LXx, Z[[m1]], qr(LXx), R, nsim = nsim, seed = seed)
  p.value <- mean(stat.sim >= stat.obs)
  
  result <- c(stat.obs = stat.obs, p.value = p.value)
  return(result[2])
}

########################################################################
###
###  Script created by Dr. Reeder, sent on Aug 15, 2016.
###  Modified as a function to run with different options
###
###  Used in the simulations, linear test comp by MYT, MAC, AMS, AM.
########################################################################

HorvathEtal<-function(Yi,Xit){
  
  Xc<-apply(Xit,2,function(x)x-mean(x))
  Y<- Yi
  N<-length(Y)
  M<-ncol(Xc)
  L<-10
  #	'basis' is a basis used to transform discrete data to functional data.
  #basis=create.fourier.basis(c(0,1),L)
  basis = create.bspline.basis(rangeval=c(0,1), nbasis=L, norder=4)
  
  #	'functional_Xc' is the functional version of 'Xc.'
  functional_Xc=Data2fd(t(Xc),argvals=seq(from=0,to=1,length.out=M),basis)
  
  d=3
  #	We now do functional PCA.
  pca=pca.fd(functional_Xc,nharm=d)
  
  #	'v' is the eigenfunctions associated with X(t).
  v=pca$harmonics
  
  r=d*(d+1)/2	#	'r' is the length of 'Ahat.'
  
  #############################################
  vech = function(A)
  {
    p=dim(A)[2]
    output = matrix(0,p*(p+1)/2,1)
    index = 0
    for (j in 1:p)
    {
      for (i in 1:j)
      {
        index = index + 1
        output[index,1] = A[i,j]
      }
    }
    output
  }
  
  vechinverse = function(x)
  {
    r=dim(x)[1]
    p=(-1 + sqrt(1+8*r) )/2
    output = matrix(0,p,p)
    rm(r)
    index = 0
    for (j in 1:p)
    {
      for (i in 1:j)
      {
        index = index + 1
        output[i,j] = x[index,1]
      }
    }
    rm(p)
    output
  }
  #########################################
  Dhat = matrix(0,r,N)
  for (n in 1:N)
  {
    Dmatrix = matrix(0,d,d)
    for (j in 1:d)
    {
      for (i in 1:j)
      {
        Dmatrix[i,j] = inprod(functional_Xc[n],v[i]) * inprod(functional_Xc[n],v[j])
      }
    }
    Dhat[,n] = vech(Dmatrix)
  }
  rm(i,j)
  
  Fhat = matrix(0,d,N)
  for (n in 1:N)
  {
    for (j in 1:d)
    {
      Fhat[j,n] = inprod(functional_Xc[n],v[j])
    }
  }
  rm(j,n)
  
  Zhat = matrix(0,N,r+d+1)
  for (n in 1:N)
  {
    Zhat[n,1:r] = t(Dhat[,n])
    Zhat[n,(r+1):(r+d)] = t(Fhat[,n])
    Zhat[n,r+d+1] = 1
  }
  rm(n)
  
  estimator = solve(t(Zhat) %*% Zhat) %*% t(Zhat) %*% Y
  uhat = estimator[r+d+1,1]
  Ahat = matrix(estimator[1:r,1],r,1)
  Ahatmatrix = vechinverse(Ahat)
  Bhat = estimator[(r+1):(r+d),1]
  Mhat = matrix(rowMeans(Dhat),r,1)
  Ghat = Dhat[,1] %*% t(Dhat[,1])
  
  for (n in 2:N)
  {
    Ghat = Ghat + Dhat[,n] %*% t(Dhat[,n])
  }
  rm(n)
  Ghat = Ghat/N
  
  sum1 = matrix(rep(0,N),N,1)
  for (n in 1:N)
  {
    for (i in 1:d)
    {
      sum1[n] = sum1[n] + Bhat[i]*inprod(functional_Xc[n],v[i])
    }
  }
  rm(i,n)
  
  sum2 = matrix(rep(0,N),N,1)
  for (n in 1:N)
  {
    for (i in 1:d)
    {
      for (j in i:d)
      {
        sum2[n] = sum2[n] + Ahatmatrix[i,j] * inprod(functional_Xc[n],v[i]) * inprod(functional_Xc[n],v[j])
      }
    }
  }
  rm(i,j,n)
  
  Ehat = Y - uhat - sum1 - sum2
  rm(sum1,sum2)
  tau_hat_squared = mean(Ehat^2)
  test_statistic = N/tau_hat_squared * t(Ahat) %*% ( Ghat - Mhat%*%t(Mhat) ) %*% Ahat
  
  results = 1-pchisq(test_statistic,df=r)
  return(results)
}

