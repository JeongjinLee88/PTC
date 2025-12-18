##  Mapping reals to positive reals
Trans = function(y){
  x=log(exp(y)+1)
  x[!is.finite(x)]=y[!is.finite(x)]
  return(x)
}
##  Mapping positive reals to reals
InvT = function(x){
  y=log(exp(x)-1)
  y[!is.finite(y) & x > 1]=x[!is.finite(y) & x > 1]
  return(y)
}
##  Transformed matrix multiplication
Amul <- function(A, x){
  y=InvT(x)
  Ay=A%*%y
  Ay=Trans(Ay)
  return(Ay)
}
##  Transformed scalar multiplication
Cmul <- function(c, x){
  y <- InvT(x)
  cy <- c*y
  cx <- Trans(cy)
  return(cx)
}
##  Transformed addition
vSum <- function(v1, v2){
  y1 <- InvT(v1)
  y2 <- InvT(v2)
  sumX <- y1 + y2
  sumV <- Trans(sumX)
  return(sumV)
}
##  Transformed subtraction
vSub <- function(v1, v2){
  y1 <- InvT(v1)
  y2 <- InvT(v2)
  sumX <- y1 - y2
  sumV <- Trans(sumX)
  return(sumV)
}
