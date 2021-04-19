################################################################################################
#
# Code to generate mesh and the Moran's Eigencomponents
#
################################################################################################
################################################################################################
library(fields) ; library(mvtnorm) ; library(classInt); library(INLA)
source(file = "../source/helperFunctions.R") # Helper file with functions
load("../samples/SpatialData.RData")
################################################################################################
#
# Generate Mesh using INLA package
#
################################################################################################
mesh <- inla.mesh.2d(gridLocation, 
                     max.edge=c(0.05,0.5),
                     cutoff = 0.025,
                     offset=c(0.1, 0.1))
# Tips for mesh construction:
# Modify max.edge and cutoff parameters in the inla.mesh.2d() function to adjust mesh density
################################################################################################
# Projector Matrix 
################################################################################################
AMat <- inla.spde.make.A(mesh, loc=gridLocation)  # model-fitting locations
AMatCV <- inla.spde.make.A(mesh, loc=CVgridLocation) # validation locations
################################################################################################
# Visualize mesh and observation locations 
################################################################################################
par(mfrow=c(1,1),mar=c(2,2,2,2))
plot(mesh,main="")
points(x=mesh$loc[,1], y=mesh$loc[,2],col="black",pch=16,cex=0.5)
points(x=gridLocation[,1], y=gridLocation[,2],col="blue",pch=16,cex=1)
points(x=CVgridLocation[,1], y=CVgridLocation[,2],col="red",pch=18,cex=1.5)
mtext("Mesh: INLA",cex=2)
legend("topright", legend=c("Model Fitting (Training)" , "Validation (Test)" , "Mesh Vertices"),
       col=c("blue","red","black"), pch=c(16,16,16))
#######################################################
# Generate Moran's Basis Functions
################################################################################################
# See Hughes and Haran (2012) or Lee and Haran (2020+) for details
DMat<-diag(apply(mesh$graph$vv,1,sum)) # Diagonal Matrix with total mumber of neighbors
WeightMat<-mesh$graph$vv # Neighborhood/Weight Matrix
PrecMat<-DMat-WeightMat # ICAR precision matrix (to be used as precision matrix of basis coefficients)
Nnew<-nrow(WeightMat) # Number of mesh vertices
OrthSpace<-diag(Nnew)-(rep(1,Nnew)%*%t(rep(1,Nnew)))/Nnew
MoransOperator<-OrthSpace%*%(WeightMat%*%OrthSpace)# Moran's Operator
MoransOperatorEig<-eigen(MoransOperator) # Moran's Basis functions (i.e. Eigencomponents of the Moran's Operator)
#######################################################
# Plot first 16 Moran's Basis functions
par(mfrow=c(4,4), mar=c(2,2,2,2))
for(i in 1:16){
  plotRF(dat=MoransOperatorEig$vectors[,i], location = mesh$loc, label = paste('Eigenvector #',i,sep=""))  
}

#######################################################
# Save File
save(mesh, AMat, AMatCV, MoransOperatorEig, file="../samples/mesh.RData")
