################################################################################################
#
# Code to fit count spatial model with spatially-varying coefficients using PICAR
#
################################################################################################
################################################################################################
library(fields) ; library(mvtnorm) ; library(classInt); library(rstan)
setwd("../..") # Set Directory
source(file = "source/helperFunctions.R") # Helper file with functions
load("samples/SpatialData.RData") # Dataset
load("samples/mesh.RData") # Mesh File
load("samples/SVC_Morans.RData") # Mesh File - Use the same as the count data
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
m_norm<-stan(file="run/stan/svc_PICAR.stan",
             data = list(N=n,
                         K=k,
                         P=p,
                         y=obsModSVC,
                         X=XMat,
                         M=M,
                         MQM=MQM),
             pars = c("beta","delta","delta_B","tau","tau_B"),
             iter=iter,chains = 1)

#Extract Pertinent Information from output
totTime<-sum(get_elapsed_time(m_norm)) # Model Fitting Time
parMat<-cbind(rstan::extract(m_norm, c('beta'))[[1]], # Extract samples for Model Parameters
              rstan::extract(m_norm, c('tau'))[[1]],
              rstan::extract(m_norm, c('tau_B'))[[1]])
deltaMat<-rstan::extract(m_norm, c('delta'))[[1]] # Extract samples for the basis coefficients.
deltaMat_B<-rstan::extract(m_norm, c('delta_B'))[[1]] # Extract samples for spatially-varying coefficients
################################################################################################
# Summarize Samples
################################################################################################
summaryMat<-list()
summaryMat[[1]]<-round(summaryFunction(parMat,
                                       totTime=totTime),3)
summaryMat[[2]]<-round(summaryFunction(deltaMat,
                                       totTime=totTime),3)
summaryMat[[3]]<-round(summaryFunction(deltaMat_B,
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
mcmcDat[[1]]<-parMat; colnames(mcmcDat[[1]])<-c("beta1","beta2","beta3","tau","tauBeta")
mcmcDat[[2]]<-deltaMat
mcmcDat[[3]]<-deltaMat_B
cvSummary<-cvFunction.svc(mcmcDat=mcmcDat,AMatCV=AMatCV,mBase=mBase,
                          XMatCV=XMatCV,obsCV=obsCVSVC)
print(cvSummary[[1]]) #MSPE
################################################################################################
# Save Data
################################################################################################
save(summaryMat,parMat,deltaMat,deltaMat_B,totTime,cvSummary,file="run/stan/svcPICARstan.RData") # Save Data
