################################################################################################
#
# Code to fit count spatial model with spatially-varying coefficients using PICAR
#
################################################################################################
################################################################################################
library(fields) ; library(mvtnorm) ; library(classInt); library(nimble)
setwd("../..") # Set Directory
source(file = "source/helperFunctions.R") # Helper file with functions
source(file = "source/nimbleSource.R") # Nimble Functions
load("samples/SpatialData.RData") # Dataset
load("samples/mesh.RData") # Mesh File
load("samples/SVC_Morans.RData") # Mesh File - Use the same as the count data example
################################################################################################
# Preliminaries for NIMBLE
################################################################################################
p<-ncol(mBase) # Number of basis functions
M<-as.matrix(AMat%*%mBase) # Projected Moran's Basis functions
Q<-diag(ncol(AMat)) # Prior precision matrix for the mesh vertices. Choices are the ICAR, CAR or an identity matrix
MQM<-t(mBase)%*%(Q%*%mBase) # Precision matrix for the basis coefficients.
niter=10000 # Iterations of MCMC Algorithm
consts   <- list(n=n,p=p)
data     <- list(Z=obsModSVC,X=XMat,M=M,MQM=MQM,mn=rep(0,p))
inits    <- list(beta1=rnorm(1),beta2=rnorm(1),beta3=rnorm(1),
                 tau=2, tauBeta=2,
                 delta=rnorm(p), bDelta=rnorm(p))
################################################################################################
# Nimble Setup
################################################################################################
Rmodel<-nimbleModel(code=count_PICAR_SVC_string, data = data,   # Build Model
                    constants=consts, inits = inits)
Cpicar <- compileNimble(Rmodel) # Compile Model
picarConf <- configureMCMC(Rmodel, print = TRUE) # Configure the model
picarConf$printSamplers() # Display samplers
picarConf$addMonitors(c("beta1", "beta2", "beta3",
                        "tau","tauBeta","delta","bDelta")) # Indicate with variables we want to sample
picarMCMC <- buildMCMC(picarConf) # Build the MCMC sampler
CpicarMCMC <- compileNimble(picarMCMC) # Compile the MCMC Sampler
################################################################################################
# Run the MCMC algorithm
################################################################################################
pt<-proc.time()
CpicarMCMC$run(niter)
ptFinal<-proc.time()-pt
samples <- as.matrix(CpicarMCMC$mvSamples)
################################################################################################
# Summarize Samples
################################################################################################
deltaInd<-grep("delta",colnames(samples)) # Indices for delta
betaDeltaInd<-grep("bDelta",colnames(samples)) # Indices for delta
summaryMat<-list()
summaryMat[[1]]<-round(summaryFunction(samples[,-c(deltaInd,betaDeltaInd)],
                                       totTime=ptFinal[3]),3)
summaryMat[[2]]<-round(summaryFunction(samples[,deltaInd],
                                       totTime=ptFinal[3]),3)
summaryMat[[3]]<-round(summaryFunction(samples[,betaDeltaInd],
                                       totTime=ptFinal[3]),3)
summaryMat[[1]] # Summary Table
apply(summaryMat[[2]],1,mean) # Mean results for the reparameterized random effects/basis coefficients
apply(summaryMat[[3]],1,mean) # Mean results for the reparameterized SVCs

################################################################################################
# Trace Plots
################################################################################################
par(mfrow=c(3,2),mar=c(2,2,2,2))  
sampInd<-floor(seq(1,nrow(samples),length.out = 1000))
for(i in 1:3){
  plot.ts(samples[sampInd,paste("beta",i,sep="")], 
          main=paste("beta",i,sep="")); 
  abline(h=beta[i],col="red",lwd=2)
  plot(density(samples[sampInd,paste("beta",i,sep="")]),
       main=paste("beta",i,sep="")); 
  abline(v=beta[i],col="red",lwd=2)
}
################################################################################################
# Compute Mean Squared Prediction Error
################################################################################################
mcmcDat<-list()
mcmcDat[[1]]<-samples[,-c(deltaInd,betaDeltaInd)]
mcmcDat[[2]]<-samples[,deltaInd]
mcmcDat[[3]]<-samples[,betaDeltaInd]
cvSummary<-cvFunction.svc(mcmcDat=mcmcDat,AMatCV=AMatCV,mBase=mBase,
                          XMatCV=XMatCV,obsCV=obsCVSVC)
print(cvSummary[[1]]) #MSPE
################################################################################################
# Save Data
################################################################################################
save(summaryMat,samples,ptFinal,cvSummary,file="run/nimble/svcPICAR.RData") # Save Data
