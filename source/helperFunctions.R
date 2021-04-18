################################################################
#
# Libraries to load 
#
################################################################
library(mvtnorm);library(fields); library(classInt)
################################################################
# Load source file for MCMC diagnosis
source("~/Dropbox/PICAR/Code/FINISHED/source/batchmeans.R")
################################################################
# Highest Posterior Density Interval
## Using Ming-Hui Chen's paper in Journal of Computational and Graphical Stats.
################################################################
hpd <- function(samp,p=0.05){
  ## to find an approximate (1-p)*100% HPD interval from a
  ## given posterior sample vector samp
  r <- length(samp)
  samp <- sort(samp)
  rang <- matrix(0,nrow=trunc(p*r),ncol=3)
  dimnames(rang) <- list(NULL,c("low","high","range"))
  for (i in 1:trunc(p*r)) {
    rang[i,1] <- samp[i]
    rang[i,2] <- samp[i+(1-p)*r]
    rang[i,3] <- rang[i,2]-rang[i,1]
  }
  hpd <- rang[order(rang[,3])[1],1:2]
  return(hpd)
}
################################################################
#
# Functions to generate Stationary and Isotropic Covariance Matrices
#
################################################################
# Exponential Covariance Function
################################################################
expCov<-function(distMat,phi){
  exp(-distMat/phi)
}
################################################################
# Gaussian Covariance Function
################################################################
sqeCov<-function(distMat,phi){
  exp(-0.5*(distMat/phi)^2)
}
################################################################
# Matern Covariance Function with smoothness nu=2.5
################################################################
matCov<-function(distMat,phi){
  (1+(sqrt(5)*(distMat/phi))+((5*distMat^2)/(3*(phi^2))))*exp(-(sqrt(5)*(distMat/phi)))
}
################################################################
# Matern Covariance Function
################################################################
Matern <- function(d, param = c(scale = 1, range = 1, smoothness = 2)) {
  scale <- param[1]
  range <- param[2]
  smoothness <- param[3]
  if (any(d < 0))
    stop("distance argument must be nonnegative")
  d <- d / range
  d[d == 0] <- 1e-10
  rootcon<-sqrt(2*smoothness)
  con <- (2^(smoothness - 1)) * gamma(smoothness)
  con <- 1 / con
  return(scale * con * ((rootcon*d)^smoothness) * besselK(rootcon*d, smoothness))
}
################################################################
# Acceptance Rate for MCMC algorithms
################################################################
accRateFunc<-function(x){
  accRate<-(length(unique(x))-1)/(length(x)-1)
  return(accRate)
}
################################################################
# Function to Summarize the MCMC output
################################################################
summaryFunction<-function(mcmcDat, # Matrix of MCMC Samples
                          totTime, # Walltime in seconds
                          bmseThresh=0.01 # Pct of MCMC mean to compare against BMSE
                          ){
  # Summary
  summaryMat<-rbind(apply(mcmcDat,2,mean),                  # Mean
                    apply(mcmcDat,2,hpd),                   # Highest Posterior Density Interval 
                    apply(mcmcDat,2,accRateFunc),           # Acceptance Rate
                    bmmat(mcmcDat)[,2],                     # Batch means standard error
                    abs(apply(mcmcDat,2,mean))*bmseThresh,  # Percent of Posterior Mean
                    apply(mcmcDat,2,ess),                   # Effective Sample Size
                    apply(mcmcDat,2,ess)/totTime)           # Effective Samplese per second
  
  rownames(summaryMat)<-c("Mean","95%CI-Low","95%CI-High",
                          "Accept","BMSE",paste(bmseThresh,"x mean"),
                          "ESS","ESS/sec")
  return(summaryMat)
}
################################################################
#
# Automated Heuristic for Rank Selection
#
################################################################
# Rank Selection for Count Data
################################################################
MLE_FindRank_Poisson<-function(XMat,   # Design Matrix (Training)
                               XMatCV, # Design Matrix (Test)
                               dimSeq, # Sequence of ranks for the Moran's Operator
                               AMat,   # Projector Matrix (Training)
                               AMatCV, # Projector Matrix (Test)
                               obsMod, # Observations (Training)
                               obsCV,  # Observations (Test)
                               MoransOperatorEig) # Moran's Eigencomponents (Test)
  {
  CVMSPE<-vector("numeric") # Container for mean squared prediction error
  betaMatList<-list() # Container for Beta Parameters
  for(i in 1:ncol(XMat)){betaMatList[[i]]<-matrix(NA,nrow=length(dimSeq),ncol=3)}
  # Loop through 
  for(jk in 1:length(dimSeq)){
    if(jk%%50==0){print(jk)}
    keepM<-1:dimSeq[jk] #Choose ranks for Moran's Basis Functions
    mBase<-(AMat%*%MoransOperatorEig$vectors[,keepM]) # Moran's Basis functions for training data
    mBaseCV<-(AMatCV%*%MoransOperatorEig$vectors[,keepM]) # Moran's Basis functions for test data
    lm1<-glm(obsMod~0+cbind(as.matrix(mBase),XMat),family = "poisson") # Fit model using GLM
    coeffs<-lm1$coefficients # Coefficients
    estMean<-coeffs[-(1:ncol(mBase))] # Estimate of Fixed Effects
    lowCI<-estMean-1.975*sqrt(diag(vcov(lm1)))[-(1:ncol(mBase))]  # 95% CI (low)
    highCI<-estMean+1.975*sqrt(diag(vcov(lm1)))[-(1:ncol(mBase))] # 95% CI (high)
    for(k in 1:length(betaMatList)){  # Save point and interval estimates for beta
      betaMatList[[k]][jk,]<-rbind(estMean,lowCI,highCI)[,k]
      }
    predCV<-exp(cbind(as.matrix(mBaseCV),XMatCV)%*%coeffs) # Prediction for test data
    CVMSPE[jk]<-mean((predCV-obsCV)^2) # Mean Squared Prediction Error
  }
  return(list(CVMSPE,betaMatList)) # Return MSPE and estimates for beta
}
################################################################
# Rank Selection for Binary Data
################################################################
MLE_FindRank_Binary<-function(XMat,   # Design Matrix (Training)
                              XMatCV, # Design Matrix (Test)
                              dimSeq, # Sequence of ranks for the Moran's Operator
                              AMat,   # Projector Matrix (Training)
                              AMatCV, # Projector Matrix (Test)
                              obsMod, # Observations (Training)
                              obsCV,  # Observations (Test)
                              MoransOperatorEig) # Moran's Eigencomponents (Test)
{
  CVMSPE<-vector("numeric") # Container for mean squared prediction error
  betaMatList<-list() # Container for Beta Parameters
  for(i in 1:ncol(XMat)){betaMatList[[i]]<-matrix(NA,nrow=length(dimSeq),ncol=3)}
  # Loop through 
  for(jk in 1:length(dimSeq)){
    if(jk%%50==0){print(jk)}
    keepM<-1:dimSeq[jk]  #Choose ranks for Moran's Basis Functions
    mBase<-(AMat%*%MoransOperatorEig$vectors[,keepM]) # Moran's Basis functions for training data
    mBaseCV<-(AMatCV%*%MoransOperatorEig$vectors[,keepM]) # Moran's Basis functions for test data
    lm1<-glm(obsMod~0+cbind(as.matrix(mBase),XMat),family = "binomial") # Fit model using GLM
    coeffs<-lm1$coefficients # Coefficients
    estMean<-coeffs[-(1:ncol(mBase))] # Estimate of Fixed Effects
    lowCI<-estMean-1.975*sqrt(diag(vcov(lm1)))[-(1:ncol(mBase))] # 95% CI (low)
    highCI<-estMean+1.975*sqrt(diag(vcov(lm1)))[-(1:ncol(mBase))] # 95% CI (high)
    for(k in 1:length(betaMatList)){# Save point and interval estimates for beta
      betaMatList[[k]][jk,]<-rbind(estMean,lowCI,highCI)[,k]
      }
    foo<-exp(cbind(as.matrix(mBaseCV),XMatCV)%*%coeffs)  
    predCV<-foo/(1+foo)
    predCV<-ifelse(predCV>0.5,1,0) # Prediction for test data
    CVMSPE[jk]<-mean((predCV-obsCV)^2) # Mean Squared Prediction Error
  }
  return(list(CVMSPE,betaMatList)) # Return MSPE and estimates for beta
}
################################################################
# Rank Selection for Gaussian Data
################################################################
MLE_FindRank_Linear<-function(XMat,   # Design Matrix (Training)
                              XMatCV, # Design Matrix (Test)
                              dimSeq, # Sequence of ranks for the Moran's Operator
                              AMat,   # Projector Matrix (Training)
                              AMatCV, # Projector Matrix (Test)
                              obsMod, # Observations (Training)
                              obsCV,  # Observations (Test)
                              MoransOperatorEig) # Moran's Eigencomponents (Test)
  {
  CVMSPE<-vector("numeric") # Container for mean squared prediction error
  betaMatList<-list() # Container for Beta Parameters
  for(i in 1:ncol(XMat)){betaMatList[[i]]<-matrix(NA,nrow=length(dimSeq),ncol=3)}
  for(jk in 1:length(dimSeq)){
    if(jk%%50==0){print(jk)}
    keepM<-1:dimSeq[jk] #Choose ranks for Moran's Basis Functions
    mBase<-(AMat%*%MoransOperatorEig$vectors[,keepM]) # Moran's Basis functions for training data
    mBaseCV<-(AMatCV%*%MoransOperatorEig$vectors[,keepM]) # Moran's Basis functions for test data
    lm1<-glm(obsMod~0+cbind(as.matrix(mBase),XMat),family = "gaussian")  
    coeffs<-lm1$coefficients # Coefficients
    estMean<-coeffs[-(1:ncol(mBase))] # Beta_Hat estimate 
    lowCI<-estMean-1.975*sqrt(diag(vcov(lm1)))[-(1:ncol(mBase))]  # 95% CI (low)
    highCI<-estMean+1.975*sqrt(diag(vcov(lm1)))[-(1:ncol(mBase))] # 95% CI (high)
    for(k in 1:length(betaMatList)){  # Save point and interval estimates for beta
      betaMatList[[k]][jk,]<-rbind(estMean,lowCI,highCI)[,k]
    }
    predCV<-cbind(as.matrix(mBaseCV),XMatCV)%*%coeffs # Prediction for test data
    CVMSPE[jk]<-mean((predCV-obsCV)^2) # Mean Squared Prediction Error
  }
  return(list(CVMSPE,betaMatList)) # Return MSPE and estimates for beta
}
################################################################
#
# Function to predict on validation dataset
#
################################################################
# Count Data
################################################################
cvFunction.pois<-function(mcmcDat, # Matrix of posterior samples 
                          AMatCV,  # Projector Matrix (test)
                          mBase,   # Moran's Basis Functions
                          XMatCV,  # Design Matrix (test)
                          obsCV,   # Test Data
                          burnin=0.5*nrow(mcmcDat[[1]])) # Burn-in
  { 
  basisMat<-AMatCV%*%mBase # Construct Basis Functions to interpolate on test locations
  cvPred<-tcrossprod(basisMat,mcmcDat[[2]][-(1:burnin),]) # Random Effects Prediction
  betaInd<-grep("beta",colnames(mcmcDat[[1]])) # index for fixed effects
  cvPreXB<-XMatCV%*%t(mcmcDat[[1]][-(1:burnin),betaInd]) # Mean prediction
  LoglambdaPred<-cvPred+cvPreXB # Linear Predictor 
  predVal<-exp(apply(LoglambdaPred,1,mean)) # Expected Value 
  cvSummary<-mean((predVal-obsCV)^2) # Mean squared prediction error
  predValuesMat<-cbind(predVal,obsCV) # Combine
  return(list(cvSummary,predValuesMat)) # MSPE and predicted values
}
################################################################
# Binary Data
################################################################
cvFunction.bin<-function(mcmcDat, # Matrix of posterior samples 
                         AMatCV,  # Projector Matrix (test)
                         mBase,   # Moran's Basis Functions
                         XMatCV,  # Design Matrix (test)
                         obsCV,   # Test Data
                         burnin=0.5*nrow(mcmcDat[[1]])) # Burn-in
{
  basisMat<-AMatCV%*%mBase # Construct Basis Functions to interpolate on test locations
  cvPred<-tcrossprod(basisMat,mcmcDat[[2]][-(1:burnin),]) # Random Effects Prediction
  betaInd<-grep("beta",colnames(mcmcDat[[1]])) # Index for fixed effects
  cvPreXB<-XMatCV%*%t(mcmcDat[[1]][-(1:burnin),betaInd]) # Mean prediction
  foo<-apply(cvPred+cvPreXB,1,mean) # Linear Predictor 
  cvPredW<-exp(foo)/(1+exp(foo)) # Expected Value 
  predVal<-ifelse(cvPredW>0.5,1,0) # Predicted value at test locations
  cvSummary<-mean((predVal-obsCV)^2)  # Mean squared prediction error
  predValuesMat<-cbind(predVal,obsCV) # Combine
  return(list(cvSummary,predValuesMat)) # MSPE and predicted values
}
################################################################
# Gaussian Data
################################################################
cvFunction.linear<-function(mcmcDat, # Matrix of posterior samples 
                            AMatCV,  # Projector Matrix (test)
                            mBase,   # Moran's Basis Functions
                            XMatCV,  # Design Matrix (test)
                            obsCV,   # Test Data
                            burnin=0.5*nrow(mcmcDat[[1]])) # Burn-in
{
  basisMat<-AMatCV%*%mBase # Construct Basis Functions to interpolate on test locations
  cvPred<-tcrossprod(basisMat,mcmcDat[[2]][-(1:burnin),]) # Random Effects Prediction
  betaInd<-grep("beta",colnames(mcmcDat[[1]])) # Index for fixed effects
  cvPreXB<-XMatCV%*%t(mcmcDat[[1]][-(1:burnin),betaInd]) # Mean prediction
  predVal<-apply(cvPred+cvPreXB,1,mean) # Linear Predictor 
  cvSummary<-mean((predVal-obsCV)^2) # Mean squared prediction error
  predValuesMat<-cbind(predVal,obsCV) # Combine
  return(list(cvSummary,predValuesMat)) # MSPE and predicted values
}
################################################################
# Count Data with Spatially-varying Coefficients
################################################################
cvFunction.svc<-function(mcmcDat, # Matrix of posterior samples 
                          AMatCV,  # Projector Matrix (test)
                          mBase,   # Moran's Basis Functions
                          XMatCV,  # Design Matrix (test)
                          obsCV,   # Test Data
                          burnin=0.5*nrow(mcmcDat[[1]])) # Burn-in
{ 
  basisMat<-AMatCV%*%mBase # Construct Basis Functions to interpolate on test locations
  cvPred<-tcrossprod(basisMat,mcmcDat[[2]][-(1:burnin),]) # Random Effects Prediction
  cvPred_svc<-tcrossprod(basisMat,mcmcDat[[3]][-(1:burnin),]) # SVC prediction
  betaInd<-grep("beta",colnames(mcmcDat[[1]])) # index for fixed effects
  cvPreXB<-XMatCV%*%t(mcmcDat[[1]][-(1:burnin),betaInd]) # Mean prediction
  LoglambdaPred<-cvPred+XMatCV[,1]*cvPred_svc+cvPreXB # Linear Predictor 
  predVal<-exp(apply(LoglambdaPred,1,mean)) # Expected Value 
  cvSummary<-mean((predVal-obsCV)^2) # Mean squared prediction error
  predValuesMat<-cbind(predVal,obsCV) # Combine
  return(list(cvSummary,predValuesMat)) # MSPE and predicted values
}
################################################################
# Figure to Plot spatial data
################################################################
plotRF<-function(dat,rangeDat=dat,label="Plot",location,length.out=10,pch=16,cex=1){
  breaks <- seq(range(rangeDat,na.rm = TRUE)[1],
                range(rangeDat,na.rm = TRUE)[2],
                length.out=length.out)
  pal <- tim.colors(length(breaks)-1,alpha = 1)
  fb <- classIntervals(dat, n = length(pal),
                       style = "fixed", fixedBreaks = breaks)
  col <- findColours(fb, pal)
  plot(x=location[,1],y=location[,2],col=col, pch=pch,cex=cex,
       main=label)
}

################################################################
#
# Ordered Categorical Data:
#
################################################################
# Ordinal Data Likelihood
################################################################
ordLlhd<-function(x,theta,linPred,log){
  cutlim<-c(0,theta)
  if(x==1){
    prob<-(1/(1+exp(linPred-cutlim[x])))
  }else if(x==4){
    prob<-1-(1/(1+exp(linPred-cutlim[x-1])))
  }else{
    prob<-(1/(1+exp(linPred-cutlim[x])))-(1/(1+exp(linPred-cutlim[x-1])))
  }
  if(log==TRUE){
    prob=log(prob)
  }
  return((prob))
}
################################################################
# Ordinal Data: Sampling Mechanism
################################################################
# Returns probability with input cutoffs and linear predictor 
ordProb<-function(theta,linPred){
  K<-length(theta)+1
  prob<-vector("numeric")
  prob[1]<-(1/(1+exp(linPred-theta[1]))) # Probability of Group 1
  prob[K]<-1-(1/(1+exp(linPred-theta[K-1])))
  for(i in 2:(K-1)){prob[i]<-(1/(1+exp(linPred-theta[i])))-(1/(1+exp(linPred-theta[i-1])))}
  return((prob))
}
# Prediction for Ordinal Data
ordinalSample<-function(x,theta){
  cutoff<-c(-Inf,theta,Inf)
  return(max(which(x>cutoff)))
}
# alpha2theta
alpha2Theta<-function(alpha){
  c(exp(alpha[1]),sum(exp(alpha)))
}

theta2Alpha<-function(theta){
  c(log(theta[1]),log(theta[2]-theta[1]))
}
################################################################
# Prediction for Validation/Test Data
################################################################
cvFunction.ordinal<-function(mcmcDat, # Matrix of posterior samples 
                         AMatCV,  # Projector Matrix (test)
                         mBase,   # Moran's Basis Functions
                         XMatCV,  # Design Matrix (test)
                         obsCV,   # Test Data
                         burnin=0.5*nrow(mcmcDat[[1]])) # Burn-in
{ 
  basisMat<-AMatCV%*%mBase # Construct Basis Functions to interpolate on test locations
  cvPred<-tcrossprod(basisMat,mcmcDat[[2]][-(1:burnin),]) # Random Effects Prediction
  betaInd<-grep("beta",colnames(mcmcDat[[1]])) # index for fixed effects
  cvPreXB<-XMatCV%*%t(mcmcDat[[1]][-(1:burnin),betaInd]) # Mean prediction
  linPred<-cvPred+cvPreXB # Linear Predictor
  linPred<-apply(linPred,1,mean) # Linear Predictor Mean
  theta<-c(0,apply(mcmcDat[[3]][-(1:burnin),], 2, mean)) # Cutoff Mean
  ncv<-length(linPred)
  predVal<- vector("numeric") 
  for(i in 1:ncv){ # Predict Category based on cutoffs
    predVal[i]<-ordinalSample(x=linPred[i],theta = theta)
  }
  cvSummary<-mean(1-(predVal==obsCV)) # Misclassification Error
  predValuesMat<-cbind(predVal,obsCV) # Combine
  return(list(cvSummary,predValuesMat)) # MSPE and predicted values
}
