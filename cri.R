library(MASS)
library(svd)

f_normalize_A_Vector<-function(x){
  centeredX<-scale(x,scale = FALSE)
  norm_of_centeredX<-sqrt(sum(centeredX^2))
  normalizedVector<-centeredX/norm_of_centeredX
  return(normalizedVector)
}

# Z: Minimal transformation + A: identity mapping
OneTimeCRI.Z = function(Xtrain, Ytrain){
  Xtrain = apply(Xtrain, 2, f_normalize_A_Vector)
  Ytrain = f_normalize_A_Vector(Ytrain)
  SVD.x = svd(Xtrain)
  cri.z = (SVD.x$v %*% t(SVD.x$u) %*% Ytrain)^2
  cri.z = as.vector(cri.z)
  names(cri.z) = colnames(Xtrain)
  return( cri.z )
}

# Z: Minimal transformation + A: correlation based reallocation 
OneTimeGCD = function(Xtrain, Ytrain){
  Xtrain = apply(Xtrain, 2, f_normalize_A_Vector)
  Ytrain = f_normalize_A_Vector(Ytrain)
  r = rankMatrix(Xtrain)
  p = dim(Xtrain)[2]
  SVD.x = svd(Xtrain, nu = r, nv = r)
  regPA.raw = (SVD.x$v %*% diag(1 / SVD.x$d[1:r]) %*% t(SVD.x$v))^2
  regPA = regPA.raw / matrix(rep(1, p), nrow = p) %*% apply(regPA.raw, 2, sum)
  GCD = regPA %*% (SVD.x$v %*% t(SVD.x$u) %*% Ytrain)^2
  GCD = as.vector(GCD)
  names(GCD) = colnames(Xtrain)
  return( GCD )
}

# Z: Minimal transformation + A: correlation based reallocation 
OneTimeCRI = function(Xtrain, Ytrain){
  Xtrain = apply(Xtrain, 2, f_normalize_A_Vector)
  Ytrain = f_normalize_A_Vector(Ytrain)
  SVD.x = svd(Xtrain)
  cri = (SVD.x$v %*% diag(SVD.x$d) %*% t(SVD.x$v))^2 %*% (SVD.x$v %*% t(SVD.x$u) %*% Ytrain)^2
  cri = as.vector(cri)
  names(cri) = colnames(Xtrain)
  return( cri )
}


