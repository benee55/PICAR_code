################################################################
#
# Nimble Functions for fitting hierarchical spatial models
#
################################################################
#
# Covariance Functions
#
################################################################
# Exponential Covariance Function
################################################################
expcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- exp(-dists[i,j]/phi)
      }
    }
    return(result)
  })
################################################################
# Matern Covariance Function with smoothness 2.5
################################################################
matcov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- (1+(sqrt(5)*(dists[i,j]/phi))+((5*dists[i,j]^2)/(3*(phi^2))))*exp(-(sqrt(5)*(dists[i,j]/phi)))
      }
    }
    return(result)
  })
################################################################
# Squared Exponential Covariance Function
################################################################
sqecov <- nimbleFunction(     
  run = function(dists = double(2), phi = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    for(i in 1:n){
      for(j in 1:n){
        result[i, j] <- exp(-(dists[i,j]/phi)^2)
      }
    }
    return(result)
  })

######################################################################
#
# Full Hierarchical Spatial models - linear, binary, and counts
#
######################################################################
# Linear Spatial Model - Full Hierarchical Spatial model
######################################################################
linear_model_string <- nimbleCode({
  # Data Model
  Z[1:n] ~ dmnorm(mean = mn[1:n], cov = fullCovMat[1:n,1:n])
  # Constant and Cov Matrix
  fullCovMat[1:n,1:n]<- sigma2*covMat[1:n,1:n]+tau2*diag(n)
  mn[1:n]<-beta1*X[,1] + beta2*X[,2] + beta3*X[,3]
  covMat[1:n,1:n]<- expcov(dists[1:n,1:n],phi) # Need to use the right Covariance function
  # Process Model
  # None needed because Ws are marginalized out 
  # Parameter Model
  # Set different priors
  sigma2   ~  dinvgamma(0.2, 0.2) 
  tau2   ~  dinvgamma(0.2, 0.2) 
  phi   ~  dunif(0,1)
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
})
######################################################################
# Binary Spatial Model - Full Hierarchical Spatial model
######################################################################
binary_model_string <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    Z[i] ~ dbinom(prob = prob[i] , size = 1)
    prob[i] <- exp(W[i]+XB[i])/(1+exp(W[i]+XB[i]))
  }
  
  # Constant and Cov Matrix
  XB[1:n]<-beta1*X[,1] + beta2*X[,2] + beta3*X[,3]
  covMat[1:n,1:n]<- expcov(dists[1:n,1:n],phi)
  fullCovMat[1:n,1:n]<- sigma2*covMat[1:n,1:n]
  
  # Process Model
  W[1:n] ~ dmnorm(mean = mn[1:n], cov = fullCovMat[1:n,1:n])
  
  # Parameter Model
  sigma2   ~  dgamma(0.2, 0.2)
  phi   ~  dunif(0,1)
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
})
######################################################################
# Count Spatial Model - Full Hierarchical Spatial model
######################################################################
count_model_string <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    lambda[i] <- exp(W[i]+XB[i])
    Z[i] ~ dpois(lambda[i])
  }
  
  # Constant and Cov Matrix
  XB[1:n]<-beta1*X[,1] + beta2*X[,2] + beta3*X[,3]
  covMat[1:n,1:n]<- expcov(dists[1:n,1:n],phi)
  fullCovMat[1:n,1:n]<- sigma2*covMat[1:n,1:n]
  
  # Process Model
  W[1:n] ~ dmnorm(mean = mn[1:n], cov = fullCovMat[1:n,1:n])
  
  # Parameter Model
  sigma2   ~  dinvgamma(0.2, 0.2)
  phi   ~  dunif(0,1)
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
})

######################################################################
#
# PICAR-based Spatial models - linear, binary, and counts
#
######################################################################
######################################################################
# PICAR for Gaussian data - Model fit using nimble
######################################################################
linear_PICAR_string <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    Z[i] ~ dnorm(mean = meanV[i] , var = sigma2)
  }
  # Mean Function
  meanV[1:n]<-XB[1:n]+W[1:n]
  XB[1:n]<-beta1*X[,1] + beta2*X[,2] + beta3*X[,3]
  W[1:n]<-M[1:n,1:p]%*%delta[1:p]
  precMat[1:p,1:p]<-tau * MQM[1:p,1:p]
  
  # Process Model
  delta[1:p] ~ dmnorm(mean = mn[1:p], prec = precMat[1:p,1:p])
  # Parameter Model
  tau   ~  dgamma(0.5,2000)
  sigma2 ~  dinvgamma(0.2, 0.2) 
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
})
######################################################################
# PICAR for Binary data - Model fit using nimble
######################################################################
binary_PICAR_string <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    prob[i] <- exp(W[i]+XB[i])/(1+exp(W[i]+XB[i]))
    Z[i] ~ dbinom(prob = prob[i] , size = 1)
  }
  
  # Constant and Cov Matrix
  XB[1:n]<-beta1*X[,1] + beta2*X[,2] + beta3*X[,3]
  W[1:n]<-M[1:n,1:p]%*%delta[1:p]
  precMat[1:p,1:p]<-tau * MQM[1:p,1:p]
  # Process Model
  delta[1:p] ~ dmnorm(mean = mn[1:p], prec = precMat[1:p,1:p])
  # Parameter Model
  tau   ~  dgamma(0.5,2000)
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
})
######################################################################
# PICAR for Count data - Model fit using nimble
######################################################################
count_PICAR_string <- nimbleCode({
  # Data Model
  for(i in 1:n){
    lambda[i] <- exp(W[i]+XB[i])
    Z[i] ~ dpois(lambda[i])
  }
  # Constant and Cov Matrix
  XB[1:n]<-beta1*X[,1] + beta2*X[,2] + beta3*X[,3]
  W[1:n]<-M[1:n,1:p]%*%delta[1:p]
  precMat[1:p,1:p]<-tau * MQM[1:p,1:p]
  # Process Model
  delta[1:p] ~ dmnorm(mean = mn[1:p], prec = precMat[1:p,1:p])
  # Parameter Model
  tau   ~  dgamma(0.5,1/2000)
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
})
######################################################################
# PICAR for Count data with Spatially-Varying Coefficients - Model fit using nimble
######################################################################
count_PICAR_SVC_string <- nimbleCode({
  # Data Model
  for(i in 1:n){
    lambda[i] <- exp(W[i]+XB[i]+B[i])
    Z[i] ~ dpois(lambda[i])
  }
  # Constant and Cov Matrix
  XB[1:n]<-beta1*X[,1] + beta2*X[,2] + beta3*X[,3]
  W[1:n]<-M[1:n,1:p]%*%delta[1:p]
  B[1:n]<-X[,1]*(M[1:n,1:p]%*%bDelta[1:p])
  precMat[1:p,1:p]<-tau * MQM[1:p,1:p]
  precMat_Beta[1:p,1:p]<-tauBeta * MQM[1:p,1:p]
  # Process Model
  delta[1:p] ~ dmnorm(mean = mn[1:p], prec = precMat[1:p,1:p])
  bDelta[1:p] ~ dmnorm(mean = mn[1:p], prec = precMat_Beta[1:p,1:p])
  # Parameter Model
  tau   ~  dgamma(0.5,1/2000)
  tauBeta   ~  dgamma(0.5,1/2000)
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
})

######################################################################
# PICAR for Ordinal data - Model fit using nimble
######################################################################
ordinal_PICAR_string <- nimbleCode({
  
  # Data Model
  for(i in 1:n){
    
    Z[i] ~ dcat(prob = prob[i,1:k])
    prob[i,1] <- (1/(1+exp(linPred[i]-theta[1]))) # Probability of Group 1
    for(j in 2:(k-1)){
      prob[i,j]<-(1/(1+exp(linPred[i]-theta[j])))-(1/(1+exp(linPred[i]-theta[j-1])))
    }
    prob[i,k] <- 1-(1/(1+exp(linPred[i]-theta[k-1])))
  }
  
  # Constant and Cov Matrix
  linPred[1:n]<-XB[1:n]+W[1:n]
  XB[1:n]<-beta1*X[,1] + beta2*X[,2] + beta3*X[,3]
  W[1:n]<-M[1:n,1:p]%*%delta[1:p]
  precMat[1:p,1:p]<-tau * MQM[1:p,1:p]
  theta[1]<-0
  theta[2]<-exp(alpha[1])
  theta[3]<-exp(alpha[1])+exp(alpha[2])
  # Process Model
  delta[1:p] ~ dmnorm(mean = mn[1:p], prec = precMat[1:p,1:p])
  # Parameter Model
  tau   ~  dgamma(0.5,2000)
  beta1 ~  dnorm(0, sd=sqrt(100))
  beta2 ~  dnorm(0, sd=sqrt(100))
  beta3 ~  dnorm(0, sd=sqrt(100))
  alpha[1] ~ dnorm(2, sd=sqrt(1))
  alpha[2] ~ dnorm(4, sd=sqrt(1))
})
