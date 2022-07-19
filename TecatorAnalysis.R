library(fda.usc)

data("tecator")

source(file="Testing_Linearity.R")
source(file="Testing_Nullity.R")

## Linearuty Tests (p-values):

Y1i<-tecator$y$Fat
Y2i<-tecator$y$Water
Y3i<-tecator$y$Protein
Xit1<-tecator$absorp.fdata

a1<-GarciaEtal(Y1i,Xit1);  a2<-GarciaEtal(Y2i,Xit1);  a3<-GarciaEtal(Y3i,Xit1)
b1<-McLeanEtal(Y1i,Xit1$data);  b2<-McLeanEtal(Y2i,Xit1$data);  b3<-McLeanEtal(Y3i,Xit1$data)
c1<-HorvathEtal(Y1i,Xit1$data);  c2<-orvathEtal(Y2i,Xit1$data);  c3<-orvathEtal(Y3i,Xit1$data)

table1<-rbind(cbind(a1,a2,a3),cbind(b1,b2,b3),cbind(c1,c2,c3))

## Nullity Tests (p-values):

a1<-GarciaEtal_null(Y1i,Xit1);  a2<-GarciaEtal_null(Y2i,Xit1);  a3<-GarciaEtal_null(Y3i,Xit1)
b1<-McLeanEtal_null(Y1i,Xit1$data);  b2<-McLeanEtal_null(Y2i,Xit1$data);  b3<-McLeanEtal_null(Y3i,Xit1$data)
c1<-HorvathEtal(Y1i,Xit1$data);  c2<-orvathEtal(Y2i,Xit1$data);  c3<-orvathEtal(Y3i,Xit1$data)

table2<-rbind(cbind(d1,d2,d3),cbind(e1,e2,e3),cbind(f1,f2,f3))

### Calculate and plot \gamma(s,t) = sum_j sum_k a_jk v_k(s)*v_j(s) for Y_1, Y_2 and Y_3

Y1i<-tecator$y$Fat
Y2i<-tecator$y$Water
Y3i<-tecator$y$Protein
Xit1<-tecator$absorp.fdata

Xc<-apply(Xit1$data,2,function(x)x-mean(x))
Y<- Y1i
N<-length(Y)
M<-ncol(Xc)
L<-10
basis = create.bspline.basis(rangeval=c(0,1), nbasis=L, norder=4)
functional_Xc=Data2fd(t(Xc),argvals=seq(from=0,to=1,length.out=M),basis)
d=3
pca=pca.fd(functional_Xc,nharm=d)
v=pca$harmonics
r=d*(d+1)/2	#	'r' is the length of 'Ahat.'

#############################################
vech = function(A){
  p=dim(A)[2]
  output = matrix(0,p*(p+1)/2,1)
  index = 0
  for (j in 1:p){
    for (i in 1:j){
      index = index + 1
      output[index,1] = A[i,j]
    }}
  output}
vechinverse = function(x){
  r=dim(x)[1]
  p=(-1 + sqrt(1+8*r) )/2
  output = matrix(0,p,p)
  rm(r)
  index = 0
  for (j in 1:p){
    for (i in 1:j){
      index = index + 1
      output[i,j] = x[index,1]
    }}
  rm(p)
  output}
Dhat = matrix(0,r,N)
for (n in 1:N){
  Dmatrix = matrix(0,d,d)
  for (j in 1:d){
    for (i in 1:j){
      Dmatrix[i,j] = inprod(functional_Xc[n],v[i]) * inprod(functional_Xc[n],v[j])
    }}
  Dhat[,n] = vech(Dmatrix)}
rm(i,j)
Fhat = matrix(0,d,N)
for (n in 1:N){
  for (j in 1:d){
    Fhat[j,n] = inprod(functional_Xc[n],v[j])
  }}
rm(j,n)

dim(Fhat);dim(Dhat)

coefa<-lm(Y~cbind(t(Fhat),t(Dhat)))$coef
a<-as.vector(coefa[5:10])
amat<-matrix(c(a[1:2],a[4],a[2:3],a[5],a[4:5],a[6]),3,3)

Gammamatrix = matrix(0,d,d)
for (j in 1:d){
  for (i in 1:j){
    Gammamatrix =amat[i,j]*(v$coefs[,i])%*%t(v$coefs[,j])
  }}

data<-data.frame(cbind(as.vector(Gammamatrix),expand.grid(1:10,1:10)))
names(data)<-c("Gammast","s","t")

library(lattice)


persp(1:10, 1:10, Gammamatrix, phi = 55, theta = 45,
      xlab = "s", ylab = "t",
      main = "Gamma(s,t)"
)


wireframe(Gammast ~ s*t, data = data,
          xlab = "s Coordinate", ylab = "t Coordinate ",
          main = "Gamma(s,t)",
          drape = TRUE,
          colorkey = TRUE,
          screen = list(z = -35, x = -60)
)
