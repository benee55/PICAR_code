
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
- Generates spatial observations of varying types - Gaussian, binary, count, ordinal, and counts with spatially-varying coefficients. The underlying spatial latent field is assumed to be a Gaussian Process with a zero-valued mean function and a covariance function from the Mat'\ern class. 
- User is able to adjust the sample sizes, model parameters, and observed/validation locations
- All samples are are saved in `/samples/SpatialData.RData`

### Step Two - Mesh Construction and Moran's Basis Functions: 
- File `/run/2_setupPICAR.R`
- Creates the mesh, Moran's eigenvector basis functions, and the projector matrices for the observed and validation locations. 
- This step corresponds to section 3.1 of Lee and Haran (2021+). 
- All components are saved in file `/samples/mesh.RData`

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

### Step Four - PICAR-based Model Fitting via MCMC:

#### Implementation A - nimble:

#### Implementation B - stan:

- Adds all functions from the script, `functions.R`, found in the `./src/` folder.

### Step Four - PICAR-based Model Fitting via MCMC: 

- Loads the raw coal-fired power plant facilities data, cleans the data, and creates a data frame with relevant covariate values. After step three, the raw facility data are stored in the `./data/` folder as `AMPD_Unit_with_Sulfur_Content_and_Regulations_with_Facility_Attributes.csv`. This section saves four output RDS files - `MonthlyUnitData.RDS`, `AnnualUnitData.RDS`, `MonthlyFacilityData.RDS`, and `AnnualFacilityData.RDS` - in the `./data/` folder. 

### Step Six: 

- Cleans all data (including environmental covariates, the SO4 response variable, and facilities data), and stores them as a single raster. This raster object is saved as `./data/central-usa-data.RDS`. This raster contains all data needed for the remaining analysis. The relevant raw data sources can be found in the subfolder, `./data/`, created in Step Two.

### Step Seven: 

- Generates posterior draws (via MCMC) from the 4 models considered in the manuscript. The samples are stored as RDS files in `./output/`. **CAUTION: THIS WILL TAKE A VERY LONG TIME**. If possible, it is recommended that the individual steps 2-5 found in `./src/so4-mcmc.R` be completed in parallel, if possible.

### Step Eight:  

- Summarizes the posterior draws with Figures, Tables, and results found in Section 4 of the manuscript. The figures will be saved as PNG files in `./output/`. **CAUTION: THIS MAY TAKE UP TO 5 HOURS** (due to the creation of Figure 4c).

### Step Nine:

- Generates the three plots found within `supp-materials.pdf`. These are saved as PNGs in `./output/`.


## Output

Upon successful completion of `main.R`, the results are saved as PNGs in the `./output/` folder. Example format includes `./output/fig1a.png`, etc. These figures will look very similar, if not identical, to those found in the manuscript. Differences can be explained by small changes induced by random draws from the posterior. However, all results should be qualitatively the same.
