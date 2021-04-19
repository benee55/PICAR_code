
# Supplemental Code for PICAR: An Efficient Extendable Approach for Fitting Hierarchical Spatial Models
Authors: Seiyon B. Lee and Murali Haran

This repository contains source code and supplemental materials from the manuscript, "PICAR: An Efficient Extendable Approach for Fitting Hierarchical Spatial Models." We provide instructions for implementing PICAR for a variety of hierarchical spatial models in nimble and stan. This repository also includes code to reproduce the main figures and tables.  

## Required Packages:
The code has been tested with R version 4.0.2, "Taking Off Again."  The following R packages must be installed before the code will run successfully:

- `fields`
- `mvtnorm`
- `classInt`
- `INLA`
- `nimble`
- `rstan`

The following programming environments must be installed prior to running the code:
- `nimble`
- `stan`

## Instructions

Before running any code, make sure the required R packages and programming environments (nimble and stan) have been installed.  Set the R working directory to the location of this README file.

### Step One - Data Generation:  
- File `/run/1_generateSamples.R`
- Generates spatial observations of varying types - Gaussian, binary, count, ordinal, and counts with spatially-varying coefficients. The underlying spatial latent field is assumed to be a Gaussian Process with a zero-valued mean function and a covariance function from the Matern class. 
- User is able to adjust the sample sizes, model parameters, and observed/validation locations
- All samples are are saved in `/samples/SpatialData.RData`
- Helper functions saved in:
  + `/source/bathmeans.R` 
  + `/source/helperFunctions.R`
<img src="/samples/data.png" width="700">

### Step Two - Mesh Construction and Moran's Basis Functions: 
- File `/run/2_setupPICAR.R`
- Creates the mesh, Moran's eigenvector basis functions, and the projector matrices for the observed and validation locations. 
- This step corresponds to section 3.1 of Lee and Haran (2021+). 
- All components are saved in file `/samples/mesh.RData`
<img src="/samples/mesh.png" width="700">
<img src="/samples/eigenvectors.png" width="700">


### Step Three - Automated Heuristic for Rank Selection:
- File `3_automatedHeuristic.R`
- Selects the appropriate rank of the Moran's eigenvector basis functions. Rank selection proceeds via the automated heuristic from Section 3,3 of Lee and Haran (2021+). 
- Rank selection proceeds separately for each data type - Gaussian, binary, count, ordinal, and counts with spatially-varying coefficients.
- The chosen rank and collection of basis functions are saved in:
  + Gaussian Data: `/samples/Linear_Morans.RData`
  + Binary Data: `/samples/Binary_Morans.RData`
  + Count Data: `/samples/Count_Morans.RData`
  + Count Data with Spatially-Varying Coefficients: `/samples/SVC_Morans.RData`
  + Ordinal Data: `/samples/Ordinal_Morans.RData`
- See code for modifications for ordinal data. 

### Step Four - PICAR-based Model Fitting via Markov Chain Monte Carlo (MCMC):

#### A. Implementation using nimble:
- Fits the PICAR-based model by drawing samples from the posterior distribution via MCMC. MCMC algorithm runs through nimble. 
- Default setting runs 10k iterations of the MCMC algorithm. 
- Provides summaries of the posterior samples (posterior mean, 95% credible intervals, acceptance rate, BMSE, number of effective samples, and effective samples per second). 
- Produces trace plots for the regression coefficients and mean-squared prediction errors for the validation sample
- Nimble source code saved as: `/source/nimbleSource.R`
- Code to run MCMC algorithms and the subsequent analyses are provided in:
  + Gaussian Data: `/run/nimble/4_runLinearPICAR.R`
  + Binary Data: `/run/nimble/4_runBinaryPICAR.R`
  + Count Data: `/run/nimble/4_runCountPICAR.R`
  + Count Data with Spatially-Varying Coefficients: `/run/nimble/4_runSVCPICAR.R`
  + Ordinal Data: `/run/nimble/4_runOrdinalPICAR.R`
- Output from model fitting and data analysis:
  + Gaussian Data: `/run/nimble/binaryPICAR.RData`
  + Binary Data: `/run/nimble/binaryPICAR.RData`
  + Count Data: `/run/nimble/countPICAR.RData`
  + Count Data with Spatially-Varying Coefficients: `/run/nimble/svcPICAR.RData`
  + Ordinal Data: `/run/nimble/ordinalPICAR.RData`

#### B. Implementation using stan:
- Similar to the nimble implemenation. However, the MCMC algorithm runs through stan. 
- Provides summaries of the posterior samples (posterior mean, 95% credible intervals, acceptance rate, BMSE, number of effective samples, and effective samples per second). 
- Produces trace plots for the regression coefficients and mean-squared prediction errors for the validation sample
- Stan source code saved in:
  + Gaussian Data: `/run/stan/linear_PICAR.stan`
  + Binary Data: `/run/stan/binary_PICAR.stan`
  + Count Data: `/run/stan/count_PICAR.stan`
  + Count Data with Spatially-Varying Coefficients: `/run/stan/svc_PICAR.stan`
  + Ordinal Data: `/run/stan/ordinal_PICAR.stan`
- Code to run MCMC algorithms and the subsequent analyses are provided in:
  + Gaussian Data: `/run/stan/4_runLinearPICAR_stan.R`
  + Binary Data: `/run/stan/4_runBinaryPICAR_stan.R`
  + Count Data: `/run/stan/4_runCountPICAR_stan.R`
  + Count Data with Spatially-Varying Coefficients: `/run/stan/4_runSVCPICAR_stan.R`
  + Ordinal Data: `/run/stan/4_runOrdinalPICAR_stan.R`
- Output from model fitting and data analysis:
   + Gaussian Data: `/run/stan/binaryPICAR.RData`
  + Binary Data: `/run/stan/binaryPICAR.RData`
  + Count Data: `/run/stan/countPICAR.RData`
  + Count Data with Spatially-Varying Coefficients: `/run/stan/svcPICAR.RData`
  + Ordinal Data: `/run/stan/ordinalPICAR.RData`

