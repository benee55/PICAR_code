################################################################################################
#
# Code to generate data: Gaussian, binary, count, oridinal, and counts with spatially-varying coefficients
#
################################################################################################
library(fields) ; library(mvtnorm) ; library(classInt)
source(file = "../source/helperFunctions.R") # Helper file with functions
################################################################################################
#Parameters
################################################################################################
set.seed(2021) # Set Seed
n=1000 ; ncv=300 # Data Generation
beta=c(0,0.5,1) ; phi=0.2 ; sigma2=1 ; tau2=0.1
################################################################################################
# Generate Locations 
################################################################################################
gridLocation<-cbind(runif(n,min = 0,max = 1),runif(n,min = 0,max = 1)) # Traning Dataset
CVgridLocation<-cbind(runif(ncv,min = 0,max = 1),runif(ncv,min = 0,max = 1)) # Test Dataset
comboLocation<-rbind(gridLocation,CVgridLocation) # Combine Locations
distMatFull<-as.matrix(rdist(comboLocation)) # Combined Distance  Matrix
modInd<-1:n # Training Data Indices
CVInd<-(n+1):nrow(distMatFull) # Test Data Indices
################################################################################################
# Covariates
################################################################################################
XMat<-cbind(runif(n,-1,1),runif(n,-1,1),runif(n,-1,1)) # Design Matrix for Training Data
XMatCV<-cbind(runif(ncv,-1,1),runif(ncv,-1,1),runif(ncv,-1,1)) # Design Matrix for Test Data
################################################################################################
# Covariance Matrix
# I used a Matern Cov function with parameters phi, sigma2, and fixed nu=2.5. 
# Nu can be adjusted by replacing matCov with expCov (exponential nu=0.5), 
# sqeCov (squared exponential nu=infinity), or the Matern function for varying smoothness parameters nu
################################################################################################
CovMat<-sigma2*matCov(distMat = distMatFull, phi = phi)
################################################################################################
# Latent Gaussian Random Field
################################################################################################
gpWFull <- as.numeric(rmvnorm(n=1,mean=rep(0,nrow(CovMat)),sigma = CovMat,method = "chol")) # Spatial random field
linPredFull<-gpWFull+rbind(XMat,XMatCV)%*%beta   # Linear Predictor
pWFullLinear<-linPredFull # Expected Value for Gaussian RV
pWFullPois<-exp(linPredFull) # Expected Value for Poisson RV
pWFullBin<-exp(linPredFull)/(1+exp(linPredFull)) # Expected Value for Bernoulli RV
################################################################################################
# Observations
################################################################################################
obsFullLinear<-rnorm(n=n+ncv , mean=linPredFull , sd=sqrt(tau2)) # Linear Observations
obsFullPois<-rpois(n=n+ncv,lambda = pWFullPois) # Count Observations
obsFullBin<-rbinom(n=n+ncv, prob = pWFullBin, size=1) # Binary Observations
################################################################################################
#
# Split into Training and Test Datasets
#
################################################################################################
################################################################################################
# Training Data
################################################################################################
gpWMod<-gpWFull[modInd]# Latent Process
pWModLinear<-pWFullLinear[modInd] # Expected Value Linear
pWModPois<-pWFullPois[modInd] # Expected Value Poisson
pWModBin<-pWFullBin[modInd] # Expected Value Binary
obsModLinear<-obsFullLinear[modInd] # Observations Linear
obsModPois<-obsFullPois[modInd] # Observations Poisson
obsModBin<-obsFullBin[modInd] # Observations Binary
################################################################################################
# Test Data
################################################################################################
gpWCV<-gpWFull[CVInd]
pWCVLinear<-pWFullLinear[CVInd]
pWCVPois<-pWFullPois[CVInd]
pWCVBin<-pWFullBin[CVInd]
obsCVLinear<-obsFullLinear[CVInd]
obsCVPois<-obsFullPois[CVInd]
obsCVBin<-obsFullBin[CVInd]
################################################################################################
# # True MSPE: MSPE if we knew the true value of the parameters and the latent spatial random effects
################################################################################################
truthCVMSPELinear<-mean((pWCVLinear-obsCVLinear)^2) # Linear
truthCVMSPEPois<-mean((pWCVPois-obsCVPois)^2) # Poisson
truthCVMSPEBin<-mean((pWCVBin-obsCVBin)^2) # Binary
################################################################################################
################################################################################################
#
# Ordinal/Ordered Categorical Data
# For K categories, we estimate K − 2 cut-points.
# To avoid identifiability issues, we fix the first cutoff point at 0 (Higgs and Hoeting 2009) and estimate the the intercept. 
#
################################################################################################
# K=4 categories
theta<-c(0,1,2) # Cutoff Points for the ordinal data. 
probMat<-t(apply(linPredFull,1,ordProb,theta=theta)) # Probabilities 
obsFullOrdinal<-apply(probMat,1,sample,x = 1:4,size=1,replace=TRUE)
obsModOrdinal<-obsFullOrdinal[modInd] # Observations Ordinal - Training
obsCVOrdinal<-obsFullOrdinal[CVInd] # Observations Ordinal - Test
# Prediction
predictionOrdinal<-vector("numeric")
linPredCV<-linPredFull[CVInd]
for(i in 1:ncv){
  predictionOrdinal[i]<-ordinalSample(x=linPredCV[i],theta = theta)
}
truthMisMatch<-mean(1-(predictionOrdinal==obsCVOrdinal))
################################################################################################
################################################################################################
#
# Spatially Varying-Coefficeints
#
################################################################################################
tMat<-rbind(c(1,0.3),c(0.3,0.2)) # T-matrix for cross-correlation (See Section 2.2)
CovMat<-kronecker(tMat,CovMat) # Full Covariance Matrix. Note: we use the same covariance matrix for the other cases
set.seed(1234)
z<-t(chol(CovMat))%*%rnorm(nrow(CovMat)) # Generate latent spatial field
indW<-1:(n+ncv) # Indices for spatial random effects
indB1<-(n+ncv+1):(2*(n+ncv)) # Indices for spatially varying coefficients
gpWFullSVC <- as.numeric(z[indW])  # Vector of Spatial random effects
betaSFullSVC <- as.numeric(z[indB1])   # Vector of Spatially varying coefficients'
linPredFullSVC<-gpWFullSVC+rbind(XMat,XMatCV)%*%beta   # Linear Predictor
pWFullSVC<-exp(linPredFullSVC+c(XMat[,1],XMatCV[,1])*betaSFullSVC) # Expected Value for Poisson RV with spatially varying coefficients
obsFullSVC<-rpois(n=n+ncv,lambda = pWFullSVC) # Full Observations
pWModSVC<-pWFullSVC[modInd] #Training  Expected Value
obsModSVC<-obsFullSVC[modInd] # Training Observations
pWCVSVC<-pWFullSVC[CVInd]   #  Test Expected Value
obsCVSVC<-obsFullSVC[CVInd] # Test Observations
truthCVMSPESVC<-mean((pWCVSVC-obsCVBin)^2) # Poisson with spatially varying coefficients 
################################################################################################
# Plots
################################################################################################
par(mfrow=c(2,3), mar=c(2,2,2,2))
plotRF(dat=obsFullLinear, location = comboLocation , label="Linear Observations")
plotRF(dat=obsFullPois, location = comboLocation , label="Count Observations")
plot(x=comboLocation[,1], y=comboLocation[,2], col = obsFullBin+1 , pch=16 , main="Binary Observations")
plotRF(dat=obsFullSVC, location = comboLocation , label="Count Observations (SVC)")
plotRF(dat=obsFullOrdinal, location = comboLocation , label="Ordinal Observations")
################################################################################################
################################################################################################
#
# Save Data
#
################################################################################################
rm(CovMat, distMatFull) # Remove Large matrices
rm(list=lsf.str()) # Remove all functions
save(list = ls(.GlobalEnv),file="../samples/SpatialData.RData") # Save Data 
