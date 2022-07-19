CreateDatat2 <- function(n=100,numT=50,tfine=seq(0,1,l=50),snre=4,sig2e=NULL,sig2FLM=4,sig2ntpx=4,
                         sig2ptpx=4,beta=c(1,1,1),nxbf=5,ntbf=5, seed=NULL,
                         extendXrange=0.01,propMissing=0,sigmax=0,
                         Xtype='McLeanEtAl',J=2, errdist, ...){
  
  nfine <- numT
  dots <- list(...)
  #browser()
  #if(Xtype=='McLeanEtAl'){
  #  J <- ifelse(is.null(dots[['J']]),500,dots[['J']])# X rougher as J increases (5 and 500 in paper)
  #}else 
  if(Xtype=='PengPaul'){
    PPxeasy<- ifelse(is.null(dots[['PPxeasy']]),FALSE,dots[['PPxeasy']])  #easy or prac datasets from fpca
  }else if(Xtype=='James'){
    Xtcors <- ifelse(is.null(dots[['Xtcors']]),.9,dots[['Xtcors']])
  }else   if(Xtype=='GoldsmithEtAl'){
    J <- ifelse(is.null(dots[['J']]),10,dots[['J']])# X rougher as J increases (5 and 500 in paper)
  }
  #  if(trueF=='hill'){
  #    if(is.null(dots[['sx']])){
  #      warning('parameters \'sx\' and \'st\' should be provided so true surface has hill shape')
  #      sx <- .3
  #      st <- .3
  #    }else{
  #      sx <- dots[['sx']]
  #      st <- dots[['st']]
  #    }
  #  }
  
  if(!is.null(propMissing) & !is.null(Ji)){
    warning('Specify only one of \'propMissing\' and \'Ji\'.  \'Ji\' will be ignored.' )
    Ji <- NULL
  }  
#   if(!is.null(sigmae) & !is.null(snre)){
#     warning('Specify only one of \'sigmae\' and \'snre\'.  \'snre\' will be ignored.' )
#     snre <- NULL
#   }
#   if(!is.null(sigmax) & !is.null(snrx)){
#     warning('Specify only one of \'sigmax\' and \'snrx\'.  \'snrx\' will be ignored.' )
#     snrx <- NULL
#   }
  
  if(!is.null(seed))
    set.seed(seed)
  
#   if(is.null(tfine))
#     tfine=seq(0,1,length=nfine) #equally spaced observation times on [0,1]  
  maxT <- max(tfine)
  mu <- numeric(nfine) # mean function is zero for all methods except YaoEtAl 
  if(Xtype=='McLeanEtAl'){
    
    #bbbeta=create.bspline.basis(rangeval=c(0,1), nbasis=20,
    #                            norder=4, breaks=NULL)
    
    # compute predictor curves and response vector for both cases of eigenvalues
    # do this nsim times and write to file
    
    efMat <- cbind( cos(pi*outer(tfine/maxT,1:J)),sin(pi*outer(tfine/maxT,1:J)) )
    ind <- rep(1:J,e=2)
    ind[seq(2,2*J,by=2)] <- ind[seq(2,2*J,by=2)]+J
    efMat <- efMat[,ind]
    #betavals=as.vector(crossprod((1:J)^(-4),cosevals))
    #truebetafd=Data2fd(tfine,betavals,bbbeta)
    #plot(truebetafd)
    xi <- sqrt(2)*cbind(sapply(1:J,function(i){rnorm(n,0,2/i)}),
                        sapply(1:J,function(i){rnorm(n,0,2/i)}) )
    lambdas <- rep(4*(1:J)^-2,e=2)
    xi <- xi[,ind]
    # Xorig <- xi%*%efMat
    
  }else if(Xtype=='James'){
    cMat <- matrix(Xtcors,nfine,nfine)
    cMat[cbind(1:nfine,1:nfine)] <- 1
    efMat <- chol(cMat)
    xi <- matrix(rnorm(n*nfine),n,nfine)
    lambdas <- 1+numeric(nfine)
  }else if(Xtype=='PengPaul'){
    if(PPxeasy){
      data(easy)
      PPdat <- easy
      #rm(easy)
    }else{
      data(prac)
      PPdat <- prac
      #rm(prac)
    }
    lambdas <- PPdat$eigenv
    bbt=create.bspline.basis(range(tfine),nbasis=10)
    tfine <- seq(min(tfine),max(tfine),l=501)
    efs <- smooth.basisPar(tfine,t(PPdat$eigenf),bbt,int2Lfd(2))$fd
    efMat <- eval.fd(tfine,efs)
    xi <- sapply(lambdas,rnorm,n=n,mean=0)
    
  }else if(Xtype=='YaoEtAl'){
    mu <- tfine+sin(tfine)
    xi <- cbind(rnorm(n,0,2),rnorm(n))
    lambdas <- c(4,1)
    efMat <- cbind(-cos(pi*tfine),sin(pi*tfine))
    #Xorig <- mu+tcrossprod(xi,efMat)
  }else if(Xtype=='GoldsmithEtAl'){
    efMat <- cbind( cos(2*pi*outer(tfine/maxT,1:J)),sin(2*pi*outer(tfine/maxT,1:J)) )
    ind <- rep(1:J,e=2)
    ind[seq(2,2*J,by=2)] <- ind[seq(2,2*J,by=2)]+J
    efMat <- efMat[,ind]
    efMat <- cbind(rep(1,nfine),tfine,efMat)
    xi <- cbind( sapply(1:J,function(i){rnorm(n,0,1/i)}), sapply(1:J,function(i){rnorm(n,0,1/i)}) )
    xi <- xi[,ind]
    xi <- cbind(5*rnorm(n,0,1),.2*rnorm(n,0,1),xi)
    lambdas <- c(25,.04,rep((1:J)^-2,e=2))
    #  Xorig=sapply(1:J,function(i){rnorm(n,0,2/i)})%*%cosevals+rnorm(1,0,5)+
    #    %*%sinevals+.2*rnorm(1)*tfine
  }else if(Xtype=='HallEtAl'){  
    mu=2*sin(2*pi*tfine/maxT)/sqrt(5) 
    xi <- sqrt(2)*rnorm(n,0,1)
    efMat <- -cos(pi*tfine/maxT)
    lambdas <- 2
  }
  Xorig <- tcrossprod(xi,efMat)
  # matplot(t(Xorig),type='l')
  Xsparse <- Xorig
  if(!is.null(propMissing)){  
    Xsparse[as.logical(rbinom(n*nfine,1,propMissing))] <- NA  
  }else{
    if(length(Ji)==1){
      for(i in 1:n) Xsparse[i,-sample(nfine,Ji)] <- NA
    }else{
      for(i in 1:n) Xsparse[i,-sample(nfine,sample(Ji,1))] <- NA
    }
  }
  
  dx <- 2
  dt <- 2
  
  DDx <- crossprod({
    DD <- diag(nxbf)
    if(dx>0) for(i in 1:dx) DD <- diff(DD)
    DD
  }) 
  DDt <- crossprod({
    DD <- diag(ntbf)
    if(dt>0) for(i in 1:dt) DD <- diff(DD)
    DD
  })
  
  
  # eigen decomp of marginal penalties:
  eDDx <- eigen(DDx)
  nullx <- (nxbf-dx+1):nxbf
  if(dx>0) eDDx$values[nullx] <- 0
  
  eDDt <- eigen(DDt)
  nullt <- (ntbf-dt+1):ntbf
  if(dt>0) eDDt$values[nullt] <- 0
  
  bbt <- create.bspline.basis(range(tfine), nbasis=ntbf)
  bbtevals <- eval.basis(tfine, bbt)
  
  
  Xrange <-range(Xorig, na.rm=TRUE)
  extXrange <- Xrange + c(-extendXrange*diff(Xrange), extendXrange*diff(Xrange))
  
  bbx <- create.bspline.basis(extXrange, nbasis=nxbf)
  L <- rep(1/numT, numT)

  UDinvmat_t <- eDDt$vectors[,1:(ntbf-2)]%*%diag(1/sqrt(eDDt$values[1:(ntbf-2)]))
  Zt <- bbtevals%*%UDinvmat_t
  UDinvmat_x <- eDDx$vectors[,1:(nxbf-2)]%*%diag(1/sqrt(eDDx$values[1:(nxbf-2)]))
  Design <- t( apply(Xorig,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-	cbind( 1+numeric(numT), xvals, xvals*tfine, #=Z_0
                 Zt*xvals, Zx,Zx*tfine, Zx[,rep(1:(nxbf-2),e=ntbf-2)]*Zt[,rep(1:(ntbf-2),t=nxbf-2)]) #=Z_p
    crossprod(L,Bxi)
  }) )  
  
  nupp <- dx*dt-1
  upind <- 1:nupp
  dims <- if(FALSE) c(nupp,ntbf-2,nxbf-2,nxbf-2,(nxbf-2)*(ntbf-2)) else c(nupp,ntbf-2,2*(nxbf-2),(nxbf-2)*(ntbf-2))
  dnames <- if(FALSE) c('para','x.Zt','Zx','Zx.t','Zx.Zt') else c('para','x.Zt','onet.Zx','Zx.Zt')
  npars <- sum(dims)
  npp <- npars-nupp
  parinds <- 1:npars
  names(parinds) <- unlist( mapply(FUN=function(y,x)rep(y,e=x),y=dnames,x=dims) )
  # browser()
  colnames(Design) <- names(parinds)

  bFLM <- rnorm(n=dims[2],sd=sqrt(sig2FLM))  
  bntpx <- rnorm(n=dims[3],sd=sqrt(sig2ntpx))  
  bptpx <- rnorm(n=dims[4],sd=sqrt(sig2ptpx))  
  
  coef <- c(beta,bFLM,bntpx,bptpx)
 # browser()
  
  Fxt <- t( apply(Xorig,1,function(xvals){
    Zx <- eval.basis(xvals,bbx)%*%UDinvmat_x
    Bxi<-  cbind( 1+numeric(numT), xvals, xvals*tfine, #=Z_0
                  Zt*xvals, Zx,Zx*tfine, Zx[,rep(1:(nxbf-2),e=ntbf-2)]*Zt[,rep(1:(ntbf-2),t=nxbf-2)]) #=Z_p
    Bxi%*%coef
  }) ) 
  
  intFvals <- drop(Design%*%coef)
  if (missing(errdist)){
    if(!is.null(sig2e) & !is.null(snre))
      warning('Both sig2e and snre specified. Ignoring sig2e.')
    if(is.null(sig2e))
      sig2e <- var(intFvals)/snre
    sigmae <- sqrt( sig2e )  
    y <- intFvals + sigmae*rnorm(n) 
  }else{
    err <- errdist(n)
    sigmae <- sd(err)
    y <- intFvals + err
  }
  
  y.oracle.ufe <- y - Design[, colnames(Design)=='x.Zt']%*%bFLM
  y.oracle <- y.oracle.ufe - Design[, 1:3]%*%beta
  if(sig2ptpx==0 & sig2ntpx==0){
    XtX <- crossprod(Design[,1:(sum(dims[1:2]))])
    diagPenalty <- c( rep(0,dims[1]),rep(sig2e/sig2FLM,dims[2]) )
  }else{
    XtX <- crossprod(Design)
    diagPenalty <- c( rep(0,dims[1]),rep(sig2e/sig2FLM,dims[2]),rep(sig2e/sig2ntpx,dims[3]),
                      rep(sig2e/sig2ptpx,dims[4]))
  }
 # df <- diag( solve(XtX+diag(diagPenalty),XtX) )
 # if(sig2ptpx==0 & sig2ntpx==0) df <- c(df,rep(0,sum(dims[3:4])))
#  csdims <- cumsum(dims)
#  dfs <- numeric(3)
    
#  for(i in 1:3)
#    dfs[i] <- sum(df[(csdims[i]+1):csdims[i+1]])
  
#  Xrng <- cbind(rep(tfine, e=n), as.vector(Xorig))[chull(cbind(rep(tfine, e=n), as.vector(Xorig))),]
  if(is.null(snre)) snre <- var(intFvals)/sig2e
    
    data=list(y = y,              # response vector
              y.oracle= y.oracle,
              y.oracle.ufe=y.oracle.ufe,
              eta = intFvals,          # N-vector of integrated F's for each sample
              Fxt = Fxt,          # evaluated true surface at grid of x and t values
              coef=coef,
              sigmae=sigmae,  # error variance
              snre=snre,
              sigmax=sigmax,
              Xorig = Xorig,           # oringal functional predictor sample data
              tfine = tfine,
  #            dfs=dfs,
              sig2FLM=sig2FLM,
              sig2ntpx=sig2ntpx,
              sig2ptpx=sig2ptpx,
              bFLM=bFLM,
              bntpx=bntpx,
              bptpx=bptpx,
              Xtype= Xtype,
              efMat=efMat,   # nfine x M matrix of eigenfunctions evaluated at tfine
              xi=xi,         # n x M matrix of p.c. scores
              mu=mu,         # nfine-vector of mean function evaluations
              lambda=lambdas,
     #         Xrng=Xrng,
              full=FALSE,
              parinds=parinds,
              bbx=bbx,
              bbt=bbt,
              UDinvmat_x=UDinvmat_x,
              UDinvmat_t=UDinvmat_t,
              nullBasis=0, #deprecated
              call = match.call())    # function call
  
  return(data)
}


CreateDatat2cc <- function(n=100,nfine=50,tfine=NULL,trueF='linear',seed=NULL,
                    sigmae=NULL,snre=10,sigmax=NULL,snrx=.01,Xtype='GoldsmithEtAl',
                    extendXrange=.05,propMissing=.1,Ji=NULL,sx=1,st=1,cons=1,
                           ConvexCombo=FALSE,lambda=1,FSIM=FALSE,FSIMcons=1, errdist, ...){
  
  # set simulations parameters for the chosen method for simulating the X's if not specified
  if(!is.null(tfine)) nfine <- length(tfine)
  dots <- list(...)
  #browser()
  if(Xtype=='McLeanEtAl'){
    J <- ifelse(is.null(dots[['J']]),500,dots[['J']])# X rougher as J increases (5 and 500 in paper)
  }else if(Xtype=='PengPaul'){
    PPxeasy<- ifelse(is.null(dots[['PPxeasy']]),FALSE,dots[['PPxeasy']])  #easy or prac datasets from fpca
  }else if(Xtype=='James'){
    Xtcors <- ifelse(is.null(dots[['Xtcors']]),.9,dots[['Xtcors']])
  }else   if(Xtype=='GoldsmithEtAl'){
    J <- ifelse(is.null(dots[['J']]),10,dots[['J']])# X rougher as J increases (5 and 500 in paper)
  }
  #  if(trueF=='hill'){
  #    if(is.null(dots[['sx']])){
  #      warning('parameters \'sx\' and \'st\' should be provided so true surface has hill shape')
  #      sx <- .3
  #      st <- .3
  #    }else{
  #      sx <- dots[['sx']]
  #      st <- dots[['st']]
  #    }
  #  }
  
  if(!is.null(propMissing) & !is.null(Ji)){
    warning('Specify only one of \'propMissing\' and \'Ji\'.  \'Ji\' will be ignored.' )
    Ji <- NULL
  }  
  if(!is.null(sigmae) & !is.null(snre)){
    warning('Specify only one of \'sigmae\' and \'snre\'.  \'snre\' will be ignored.' )
    snre <- NULL
  }
  if(!is.null(sigmax) & !is.null(snrx)){
    warning('Specify only one of \'sigmax\' and \'snrx\'.  \'snrx\' will be ignored.' )
    snrx <- NULL
  }
  
  if(!is.null(seed))
    set.seed(seed)
  
  if(is.null(tfine))
    tfine=seq(0,1,length=nfine) #equally spaced observation times on [0,1]  
  maxT <- max(tfine)
  mu <- numeric(nfine) # mean function is zero for all methods except YaoEtAl 
  if(Xtype=='McLeanEtAl'){
    
    #bbbeta=create.bspline.basis(rangeval=c(0,1), nbasis=20,
    #                            norder=4, breaks=NULL)
    
    # compute predictor curves and response vector for both cases of eigenvalues
    # do this nsim times and write to file
    
    efMat <- cbind( cos(pi*outer(tfine/maxT,1:J)),sin(pi*outer(tfine/maxT,1:J)) )
    ind <- rep(1:J,e=2)
    ind[seq(2,2*J,by=2)] <- ind[seq(2,2*J,by=2)]+J
    efMat <- efMat[,ind]
    #betavals=as.vector(crossprod((1:J)^(-4),cosevals))
    #truebetafd=Data2fd(tfine,betavals,bbbeta)
    #plot(truebetafd)
    xi <- sqrt(2)*cbind(sapply(1:J,function(i){rnorm(n,0,2/i)}),
                        sapply(1:J,function(i){rnorm(n,0,2/i)}) )
    lambdas <- rep(4*(1:J)^-2,e=2)
    xi <- xi[,ind]
    # Xorig <- xi%*%efMat
    
  }else if(Xtype=='James'){
    cMat <- matrix(Xtcors,nfine,nfine)
    cMat[cbind(1:nfine,1:nfine)] <- 1
    efMat <- chol(cMat)
    xi <- matrix(rnorm(n*nfine),n,nfine)
    lambdas <- 1+numeric(nfine)
  }else if(Xtype=='PengPaul'){
    if(PPxeasy){
      data(easy)
      PPdat <- easy
      #rm(easy)
    }else{
      data(prac)
      PPdat <- prac
      #rm(prac)
    }
    lambdas <- PPdat$eigenv
    bbt=create.bspline.basis(range(tfine),nbasis=10)
    tgrid <- seq(min(tfine),max(tfine),l=501)
    efs <- smooth.basisPar(tgrid,t(PPdat$eigenf),bbt,int2Lfd(2))$fd
    efMat <- eval.fd(tfine,efs)
    xi <- sapply(lambdas,rnorm,n=n,mean=0)
    
  }else if(Xtype=='YaoEtAl'){
    mu <- tfine+sin(tfine)
    xi <- cbind(rnorm(n,0,2),rnorm(n))
    lambdas <- c(4,1)
    efMat <- cbind(-cos(pi*tfine),sin(pi*tfine))
    #Xorig <- mu+tcrossprod(xi,efMat)
  }else if(Xtype=='GoldsmithEtAl'){
    efMat <- cbind( cos(2*pi*outer(tfine/maxT,1:J)),sin(2*pi*outer(tfine/maxT,1:J)) )
    ind <- rep(1:J,e=2)
    ind[seq(2,2*J,by=2)] <- ind[seq(2,2*J,by=2)]+J
    efMat <- efMat[,ind]
    efMat <- cbind(rep(1,nfine),tfine,efMat)
    xi <- cbind( sapply(1:J,function(i){rnorm(n,0,1/i)}), sapply(1:J,function(i){rnorm(n,0,1/i)}) )
    xi <- xi[,ind]
    xi <- cbind(5*rnorm(n,0,1),.2*rnorm(n,0,1),xi)
    lambdas <- c(25,.04,rep((1:J)^-2,e=2))
    #  Xorig=sapply(1:J,function(i){rnorm(n,0,2/i)})%*%cosevals+rnorm(1,0,5)+
    #    %*%sinevals+.2*rnorm(1)*tfine
  }else if(Xtype=='HallEtAl'){  
    mu=2*sin(2*pi*tfine/maxT)/sqrt(5) 
    xi <- sqrt(2)*rnorm(n,0,1)
    efMat <- -cos(pi*tfine/maxT)
    lambdas <- 2
  }
  Xorig <- tcrossprod(xi,efMat)
  # matplot(t(Xorig),type='l')
  Xsparse <- Xorig
  if(!is.null(propMissing)){  
    Xsparse[as.logical(rbinom(n*nfine,1,propMissing))] <- NA  
  }else{
    if(length(Ji)==1){
      for(i in 1:n) Xsparse[i,-sample(nfine,Ji)] <- NA
    }else{
      for(i in 1:n) Xsparse[i,-sample(nfine,sample(Ji,1))] <- NA
    }
  }
  
  if(is.null(sigmax))
    sigmax <- sqrt(var(as.vector(Xorig))/snrx)
  errs <- sigmax*matrix(rnorm(nfine*n),n,nfine)
  Xwerr <- Xorig + errs
  Xsparse <- Xsparse + errs
  if(!ConvexCombo){
    Fxt <- t( sapply(1:n,function(i){
      GenSurface(tval=tfine,type=trueF,xval=Xorig[i,],sx=sx,st=st,cons=cons)
    }) )     
  }else{
    Fxt <- t( sapply(1:n,function(i){
      lambda*GenSurface(tval=tfine,type=trueF[1],xval=Xorig[i,],sx=sx[1],st=st[1],cons=cons[1])+
      (1-lambda)*GenSurface(tval=tfine,type=trueF[2],xval=Xorig[i,],sx=sx[2],st=st[2],cons=cons[2])
    }) )   
  }
  stopifnot(nrow(Fxt)==n)
  intFvals=rowSums(Fxt)/nfine
  
  if(FSIM=='sin'){
    intFvals <- FSIMcons*sin(intFvals)
  }else if(FSIM=='square'){
    intFvals <- FSIMcons*intFvals^2
  }else if (FSIM=='root'){
   # return(list(FSIMcons,intFvals))
    intFvals <- FSIMcons*intFvals^.5
  }
  
  if (missing(errdist)){
    if(is.null(sigmae))
      sigmae <- sqrt( var(intFvals)/snre )  
    y <- intFvals + sigmae*rnorm(n) 
  }else{
    err <- errdist(n)
    sigmae <- sd(err)
    y <- intFvals + err
  }
  
  Xrng <- cbind(rep(tfine, e=n), as.vector(Xorig))[chull(cbind(rep(tfine, e=n), as.vector(Xorig))),]
  
  data=list(y = y,              # response vector
            eta = intFvals,          # N-vector of integrated F's for each sample
            Fxt = Fxt,          # evaluated true surface at grid of x and t values
            sigmae=sigmae,  # sqrt of error variance
            sigmax=sigmax,
            Xorig = Xorig,           # oringal functional predictor sample data
            tgrid = tfine,
            Xsparse = Xsparse,   # sparsified functional predictors
            Xwerr = Xwerr,
            Xtype= Xtype,
            efMat=efMat,   # nfine x M matrix of eigenfunctions evaluated at tgrid
            xi=xi,         # n x M matrix of p.c. scores
            mu=mu,         # nfine-vector of mean function evaluations
            lambda=lambdas,
            Xrng=Xrng,
            trueF = trueF,    # true surface F(x,t)
            sx = sx,
            st = st,
            cons=cons,
            call = match.call())    # function call
  
  return(data)
}

GenSurface=function(tval=seq(0,1,l=100),type,xval=NULL,xfd=NULL,mx=0,mt=.5,sx=5,st=.3,cons=1){
  if(length(xval)==0)
    xval=eval.fd(tval,xfd)
  if(type=='linear'){
    #return(cons*eval.fd(tval,truebetafd)*xval)
    return(cons*tval*xval)
  }else if (type=='square'){
    return( cons*(tval*xval)^2 )
  }else if (type=='linear2'){
    return(cons*2*sin(pi*tval)*xval)
  }else if(type=='sinlin'){
    return(sin(cons*tval*xval/sx))
  }else if(type=='cos'){
    return(cons*cos(-xval/sx+tval/st-1))
  }else if(type=='cos2'){
    return(cons*2*cos(2-xval/sx-tval/st))
  }else if(type=='sin'){
    return(cons*2*sin(2-xval/sx-tval/st))
  }else if(type=='exp'){
    return(cons*tval/st*exp(xval/sx))
  }else if(type=='hill'){
    return( cons*( exp(-(xval-mx)^2/sx^2 - (tval - 
                                              mt)^2/st^2)-.5 ))
  }else if(type=='Wood'){
    return(
      cons*(pi*sx * st) * (1.2 * exp(-(xval - mx)^2/sx^2 - (tval - 
                                                              mt)^2/st^2) + 0.8 * exp(-(xval +2)^2/sx^2 - 
                                                                                        (tval - mt)^2/st^2))
    )
  }
}

PlotData <- function(data){
  with(data, {
  x.values <- seq(min(Xorig), max(Xorig), l=25)
  t.values <- seq(min(tgrid), max(tgrid), l=25)
  grid.vals <- expand.grid(t.values, x.values)
  x <- grid.vals$Var2
  y <- grid.vals$Var1
  lam <- eval(data$call$lambda)
  
  Fxt <- lam*GenSurface(tval=y, type=trueF[1], xval=x, sx=sx[1],
                           st=st[1], cons=cons[1])+
      (1-lam)*GenSurface(tval=y, type=trueF[2], xval=x, sx=sx[2], st=st[2],
                            cons=cons[2])
  wireframe(Fxt~x*y, xlab='x', ylab='t', zlab='F(x,t)')     
  })  # end with
  
}


CreateDataFSIM <- function(n=100,nfine=50,tfine=NULL,trueF='linear',seed=NULL,
                           r.parm=.1, FSIM=FALSE,...){
  
  # set simulations parameters for the chosen method for simulating the X's if not specified
  if(!is.null(tfine)) nfine <- length(tfine)
  dots <- list(...)
  #browser()

  #  if(trueF=='hill'){
  #    if(is.null(dots[['sx']])){
  #      warning('parameters \'sx\' and \'st\' should be provided so true surface has hill shape')
  #      sx <- .3
  #      st <- .3
  #    }else{
  #      sx <- dots[['sx']]
  #      st <- dots[['st']]
  #    }
  #  }

  
  if(!is.null(seed))
    set.seed(seed)
  
  if(is.null(tfine))
    tfine=seq(0,1,length=nfine) #equally spaced observation times on [0,1]  
  maxT <- max(tfine)
  mu <- tfine # mean function is zero for all methods except YaoEtAl 
    
    #bbbeta=create.bspline.basis(rangeval=c(0,1), nbasis=20,
    #                            norder=4, breaks=NULL)
    
    # compute predictor curves and response vector for both cases of eigenvalues
    # do this nsim times and write to file
 # browser()
    efMat <- cbind( cos(pi*outer(tfine/maxT,c(2, 4))),sin(pi*outer(tfine/maxT,c(2, 4))) )
    ind <- rep(1:2,e=2)
    ind[seq(2,2*2,by=2)] <- ind[seq(2,2*2,by=2)]+2
    efMat <- efMat[,ind]
    #betavals=as.vector(crossprod((1:J)^(-4),cosevals))
    #truebetafd=Data2fd(tfine,betavals,bbbeta)
    #plot(truebetafd)
  lambdas <- c(0, .5, .25, .125)
    xi <- 1/sqrt(2)*cbind(sapply(lambdas,function(i){rnorm(n,0,i)}),
                        sapply(lambdas,function(i){rnorm(n,0,i)}) )

    xi <- xi[,ind]
    # Xorig <- xi%*%efMat

 
  Xorig <- tcrossprod(xi,efMat)+mu
  # matplot(t(Xorig),type='l')

  betavals <- 1/sqrt(3)*efMat[, 1] + 1/sqrt(3)*efMat[, 1] + 1/sqrt(6)*efMat[, 3] + 1/sqrt(6)*efMat[, 4]
  Fxt <- apply(Xorig, 2, function(x)x*betavals)
  
  stopifnot(nrow(Fxt)==n)
  intFvals=rowSums(Fxt)/nfine
  
  if(FSIM=='cos'){
    intFvals <- cos(intFvals)
  }else if(FSIM=='square'){
    intFvals <- intFvals^2
  }else if (FSIM=='root'){
    # return(list(FSIMcons,intFvals))
    intFvals <- intFvals^.5
  }else if(FSIM=='sin'){
    intFvals <- sin(intFvals)    
  }else if(FSIM){
    stop('Invalid value for FSIM specified')
  }
  
  sigmae <- sqrt(r.parm*var(intFvals))
  y <- intFvals + rnorm(n) 
  
  Xrng <- cbind(rep(tfine, e=n), as.vector(Xorig))[chull(cbind(rep(tfine, e=n), as.vector(Xorig))),]
  
  data=list(y = y,              # response vector
            eta = intFvals,          # N-vector of integrated F's for each sample
            Fxt = Fxt,          # evaluated true surface at grid of x and t values
            sigmae=sigmae,  # sqrt of error variance
            Xorig = Xorig,           # oringal functional predictor sample data
            tgrid = tfine,
            efMat=efMat,   # nfine x M matrix of eigenfunctions evaluated at tgrid
            xi=xi,         # n x M matrix of p.c. scores
            mu=mu,         # nfine-vector of mean function evaluations
            lambda=lambdas,
            Xrng=Xrng,
            betavals=betavals,
            call = match.call())    # function call
  
  return(data)
}