################################################################################################
#
# Code to fit count spatial model using PICAR
#
################################################################################################
################################################################################################
library(fields) ; library(mvtnorm) ; library(classInt); library(rstan)
setwd("../..") # Set Directory
source(file = "source/helperFunctions.R") # Helper file with functions
load("samples/SpatialData.RData") # Dataset
load("samples/mesh.RData") # Mesh File
load("samples/Linear_Morans.RData") # Mesh File
################################################################################################
# Preliminaries for NIMBLE
################################################################################################
# Inputs for Stan
k=ncol(XMat) # Design Matrix dimensions 
p<-ncol(mBase) # Number of basis functions
M<-as.matrix(AMat%*%mBase) # Projected Moran's Basis functions
Q<-diag(ncol(AMat)) # Prior precision matrix for the mesh vertices. Choices are the ICAR, CAR or an identity matrix
MQM<-t(mBase)%*%(Q%*%mBase) # Precision matrix for the basis coefficients.
################################################################################################
# Stan Setup
################################################################################################
# Fit Model Using Stan
iter=10000
m_norm<-stan(file="run/stan/linear_PICAR.stan",
             data = list(N=n,
                         K=k,
                         P=p,
                         y=obsModLinear,
                         X=XMat,
                         M=M,
                         MQM=MQM),
             pars = c("beta","delta","tau","sigma2"),
             iter=iter,chains = 1)

#Extract Pertinent Information from output
totTime<-sum(get_elapsed_time(m_norm)) # Model Fitting Time
parMat<-cbind(rstan::extract(m_norm, c('beta'))[[1]], # Extract samples for Model Parameters
              rstan::extract(m_norm, c('tau'))[[1]],
              rstan::extract(m_norm, c('sigma2'))[[1]])
deltaMat<-rstan::extract(m_norm, c('delta'))[[1]] # Extract samples for the basis coefficients.
################################################################################################
# Summarize Samples
################################################################################################
summaryMat<-list()
summaryMat[[1]]<-round(summaryFunction(parMat,
                                       totTime=totTime),3)
summaryMat[[2]]<-round(summaryFunction(deltaMat,
                                       totTime=totTime),3)
summaryMat[[1]] # Summary Table
apply(summaryMat[[2]],1,mean) # Mean results for the basis coefficients.

################################################################################################
# Trace Plots
################################################################################################
par(mfrow=c(3,2),mar=c(2,2,2,2))  
sampInd<-floor(seq(1,nrow(parMat),length.out = 1000))
for(i in 1:3){
  plot.ts(parMat[sampInd,i], 
          main=paste("beta",i,sep="")); 
  abline(h=beta[i],col="red",lwd=2)
  plot(density(parMat[sampInd,i]),
       main=paste("beta",i,sep="")); 
  abline(v=beta[i],col="red",lwd=2)
}
################################################################################################
# Compute Mean Squared Prediction Error
################################################################################################
mcmcDat<-list()
mcmcDat[[1]]<-parMat; colnames(mcmcDat[[1]])<-c("beta1","beta2","beta3","tau","sigma2")
mcmcDat[[2]]<-deltaMat
cvSummary<-cvFunction.linear(mcmcDat=mcmcDat,AMatCV=AMatCV,mBase=mBase,
                           XMatCV=XMatCV,obsCV=obsCVLinear)
print(cvSummary[[1]]) #MSPE
################################################################################################
# Save Data
################################################################################################
save(summaryMat,parMat,deltaMat,totTime,cvSummary,file="run/stan/linearPICARstan.RData") # Save Data
