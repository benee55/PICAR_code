################################################################################################
#
# Code to perform the automated rank selection
#
################################################################################################
################################################################################################
library(fields) ; library(mvtnorm) ; library(classInt); library(INLA)
source(file = "../source/helperFunctions.R") # Helper file with functions
load("../samples/SpatialData.RData") # Dataset
load("../samples/mesh.RData") # Mesh File
################################################################################################
#
# Automated heuristic for rank selection. 
# We demonstrate on Gaussian, binary, count, and count with spatially varying coefficients
#
################################################################################################
#
# 1. Gaussian Data
#
################################################################################################
# Select Rank for PICAR
################################################################################################
dimSeq<-seq(2,300,by=10) # Choose the range for the rank selection search
heuristicResults<-MLE_FindRank_Linear(XMat=XMat,XMatCV=XMatCV,dimSeq=dimSeq, # Computes MSPE for each rank chosen
                                      AMat=AMat,AMatCV=AMatCV,obsCV=obsCVLinear,obsMod=obsModLinear,
                                      MoransOperatorEig = MoransOperatorEig)
################################################################################################
# Visualize Results
################################################################################################
par(mfrow=c(2,2))
plot(x=dimSeq,y=heuristicResults[[1]],typ="l" , main="MSPE") # MSPE vs. rank
for(k in 1:ncol(XMat)){
  plot(x=dimSeq,y=heuristicResults[[2]][[k]][,1],typ="l",ylim=range(heuristicResults[[2]][[k]],na.rm = TRUE),
       main=paste("Beta",k))  # B1 estimate and confidence intervals vs. rank
  lines(x=dimSeq,y=heuristicResults[[2]][[k]][,2],lty=2) # B2 estimate and confidence intervals vs. rank
  lines(x=dimSeq,y=heuristicResults[[2]][[k]][,3],lty=2) # B3 estimate and confidence intervals vs. rank
}
################################################################################################
# Select Rank of Moran's Basis Functions. 
# Choose rank that yields the lowest MSPE
pBase<-dimSeq[which.min(heuristicResults[[1]])] ; print(pBase)
mBase<-MoransOperatorEig$vectors[,1:pBase] # Moran's Basis Function
save(pBase,mBase, file="../samples/Linear_Morans.RData")
################################################################################################
#
# 2. Binary Data
#
################################################################################################
#Select Rank for PICAR
################################################################################################
dimSeq<-seq(2,100,by=1) # Choose the range for the rank selection search
heuristicResults<-MLE_FindRank_Binary(XMat=XMat,XMatCV=XMatCV,dimSeq=dimSeq, # Computes MSPE for each rank chosen
                                      AMat=AMat,AMatCV=AMatCV,obsCV=obsCVBin,obsMod=obsModBin,
                                      MoransOperatorEig = MoransOperatorEig)
################################################################################################
# Visualize Results
################################################################################################
par(mfrow=c(2,2))
plot(x=dimSeq,y=heuristicResults[[1]],typ="l" , main="MSPE") # MSPE vs. rank
for(k in 1:ncol(XMat)){
  plot(x=dimSeq,y=heuristicResults[[2]][[k]][,1],typ="l",ylim=range(heuristicResults[[2]][[k]],na.rm = TRUE),
       main=paste("Beta",k))  # B1 estimate and confidence intervals vs. rank
  lines(x=dimSeq,y=heuristicResults[[2]][[k]][,2],lty=2) # B2 estimate and confidence intervals vs. rank
  lines(x=dimSeq,y=heuristicResults[[2]][[k]][,3],lty=2) # B3 estimate and confidence intervals vs. rank
}
################################################################################################
# Select Rank of Moran's Basis Functions. 
# Choose rank that yields the lowest MSPE
pBase<-dimSeq[which.min(heuristicResults[[1]])] ; print(pBase)
mBase<-MoransOperatorEig$vectors[,1:pBase] # Moran's Basis Function
save(pBase,mBase, file="../samples/Binary_Morans.RData")
################################################################################################
#
# 3. Count Data 
#
################################################################################################
#Select Rank for PICAR
################################################################################################
dimSeq<-seq(2,100,by=1) # Choose the range for the rank selection search
heuristicResults<-MLE_FindRank_Poisson(XMat=XMat,XMatCV=XMatCV,dimSeq=dimSeq, # Computes MSPE for each rank chosen
                                      AMat=AMat,AMatCV=AMatCV,obsCV=obsCVPois,obsMod=obsModPois,
                                      MoransOperatorEig = MoransOperatorEig)
################################################################################################
# Visualize Results
################################################################################################
par(mfrow=c(2,2))
plot(x=dimSeq,y=heuristicResults[[1]],typ="l" , main="MSPE") # MSPE vs. rank
for(k in 1:ncol(XMat)){
  plot(x=dimSeq,y=heuristicResults[[2]][[k]][,1],typ="l",ylim=range(heuristicResults[[2]][[k]],na.rm = TRUE),
       main=paste("Beta",k))  # B1 estimate and confidence intervals vs. rank
  lines(x=dimSeq,y=heuristicResults[[2]][[k]][,2],lty=2) # B2 estimate and confidence intervals vs. rank
  lines(x=dimSeq,y=heuristicResults[[2]][[k]][,3],lty=2) # B3 estimate and confidence intervals vs. rank
}
################################################################################################
# Select Rank of Moran's Basis Functions. 
# Choose rank that yields the lowest MSPE
pBase<-dimSeq[which.min(heuristicResults[[1]])] ; print(pBase)
mBase<-MoransOperatorEig$vectors[,1:pBase] # Moran's Basis Function
save(pBase,mBase, file="../samples/Count_Morans.RData")
################################################################################################
#
# 4. Count Data with Spatially-Varying Coefficients
# Rank selection commences similar to the Count case
#
################################################################################################
#Select Rank for PICAR
################################################################################################
dimSeq<-seq(2,100,by=1) # Choose the range for the rank selection search
heuristicResults<-MLE_FindRank_Poisson(XMat=XMat,XMatCV=XMatCV,dimSeq=dimSeq, # Computes MSPE for each rank chosen
                                       AMat=AMat,AMatCV=AMatCV,obsCV=obsCVSVC,obsMod=obsModSVC,
                                       MoransOperatorEig = MoransOperatorEig)
################################################################################################
# Visualize Results
################################################################################################
par(mfrow=c(2,2))
plot(x=dimSeq,y=heuristicResults[[1]],typ="l" , main="MSPE") # MSPE vs. rank
for(k in 1:ncol(XMat)){
  plot(x=dimSeq,y=heuristicResults[[2]][[k]][,1],typ="l",ylim=range(heuristicResults[[2]][[k]],na.rm = TRUE),
       main=paste("Beta",k))  # B1 estimate and confidence intervals vs. rank
  lines(x=dimSeq,y=heuristicResults[[2]][[k]][,2],lty=2) # B2 estimate and confidence intervals vs. rank
  lines(x=dimSeq,y=heuristicResults[[2]][[k]][,3],lty=2) # B3 estimate and confidence intervals vs. rank
}
################################################################################################
# Select Rank of Moran's Basis Functions. 
# Choose rank that yields the lowest MSPE
pBase<-dimSeq[which.min(heuristicResults[[1]])] ; print(pBase)
mBase<-MoransOperatorEig$vectors[,1:pBase] # Moran's Basis Function
save(pBase,mBase, file="../samples/SVC_Morans.RData")
################################################################################################
#
# 5. Ordered Categorical (Ordinal) Data
#
################################################################################################
#Select Rank for PICAR
# We treat this as a binary problem. We combine categories 1 and 2 as well as 3 and 4. 
################################################################################################
newObsMod<-ifelse(obsModOrdinal<=2,0,1) 
newObsCV<-ifelse(obsCVOrdinal<=2,0,1)
dimSeq<-seq(2,100,by=1) # Choose the range for the rank selection search
heuristicResults<-MLE_FindRank_Binary(XMat=XMat,XMatCV=XMatCV,dimSeq=dimSeq, # Computes MSPE for each rank chosen
                                      AMat=AMat,AMatCV=AMatCV,obsCV=newObsCV,obsMod=newObsMod,
                                      MoransOperatorEig = MoransOperatorEig)
################################################################################################
# Visualize Results
################################################################################################
par(mfrow=c(2,2))
plot(x=dimSeq,y=heuristicResults[[1]],typ="l" , main="MSPE") # MSPE vs. rank
for(k in 1:ncol(XMat)){
  plot(x=dimSeq,y=heuristicResults[[2]][[k]][,1],typ="l",ylim=range(heuristicResults[[2]][[k]],na.rm = TRUE),
       main=paste("Beta",k))  # B1 estimate and confidence intervals vs. rank
  lines(x=dimSeq,y=heuristicResults[[2]][[k]][,2],lty=2) # B2 estimate and confidence intervals vs. rank
  lines(x=dimSeq,y=heuristicResults[[2]][[k]][,3],lty=2) # B3 estimate and confidence intervals vs. rank
}
################################################################################################
# Select Rank of Moran's Basis Functions. 
# Choose rank that yields the lowest MSPE
pBase<-dimSeq[which.min(heuristicResults[[1]])] ; print(pBase)
mBase<-MoransOperatorEig$vectors[,1:pBase] # Moran's Basis Function
save(pBase,mBase, file="../samples/Ordinal_Morans.RData")
