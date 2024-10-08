
# Computational material for "Visibility Graph-Based Covariance Functions for Scalable Spatial Analysis in Non-Convex Partially Euclidean Domains"
 (https://arxiv.org/abs/2307.11941)

## Overview

This repository contains the code used in the paper titled "Visibility graph-based covariance functions for scalable spatial analysis in non-convex partially Euclidean domains" by Brian Gilbert and Abhirup Datta.

Data for the Chesapeake Bay data analysis is available online courtesy of the Chesapeake Bay Project (see main text for citation) and can be obtained by contacting the first author.

### Notes

- `functions.R` in the main directory implements the visGP method.
- Many '.R' files are designed to be run in a batch scripting environment. For the simulation studies, this is often to pass in a random seed. To reproduce the entire study, the scripts would be run for the seeds 1 through 4500.
- As of the time of writing, the boraGP package is available at https://github.com/jinbora0720/boraGP. We use a minor adjustment to avoid problems when barriers have zero width, which occurs in the U-shaped domain. This minor adjustment is implemented at https://github.com/bjg345/boraGP.

### Directories

#### chesapeake 

Reproduces the Chesapeake Bay data analysis.

- **chesapeake.R**: Conducts a 3-fold cross-validation analysis.
- **chesapeake_loo.R**: Performs leave-one-out cross-validation.
- **summarize.R**: Generates figures and tables from the analysis outputs.
#### sim_data_gen

This directory contains scripts for generating simulation data for the fork-shaped and U-shaped domains.

The files are organized by domain shape (`fork` or `U`) and should be run in the following sequence:

1. **points** scripts (i.e., `gen_fork_points.R`, `gen_U_points.R`) - Generate points within the specified domains.
2. **A** scripts (i.e., `gen_fork_A.R`, `gen_U_A.R`) - Create adjacency matrices for the points.
3. **ground** scripts (i.e., `gen_fork_ground.R`, `gen_U_ground.R`) - Generate "ground truth" (mean) function values.
4. **data** scripts (i.e., `gen_fork_data.R`, `gen_U_data.R`) - Create the .Rdata files used for analysis.

- **write_grid.R**: Converts the `.rds` files generated by the above scripts into `.csv` format for use with the GLGP method, which is implemented in MATLAB.
- **make_graphs.R**: Plots the fork-shaped domain.

#### fork_sim 
The `fork_sim` directory contains scripts for running various simulation methods on the fork-shaped domain data. This includes several predictive modeling methods and experiments detailed in the study's main text and supplement.
- **bora.R**: Runs the boraGP method.
- **visgp_fit.R**: Fits the visGP model.
- **maxprec.R**, **nearest.R**, **weighted.R**, **euclid.R**: Implements different visGP prediction methods (maxprec is detailed in the main paper, others in the supplement). "maxprec" refers to "maximum precision," "nearest" to "nearest clique," "weighted" to "precision-weighted," and "euclid" to standard kriging.
- **BRISC.R**: Fits the BRISC model.
- **fit.m**: Fits the GLGP model using MATLAB (dependent on other `.m` files in the directory, provided by creators of GLGP; see main text for citation).
- **gather.R**: Collects results from all simulation runs.
- **simplot.R**: Generates plots from the simulation results.

#### Subdirectories

- **dense**: Contains files analogous to those described above pertaining to the dense holdout set experiment described in the supporting information.
- **addedge**: Includes `visgp_fit_addedge.R` and `addedge_summary.R` for the experiment in the supplement investigating the effects of chordal completion on the adjacency graph.
#### u_sim

The `u_sim` directory focuses on simulations for the U-shaped domain, specifically designed for parameter and likelihood comparisons between boraGP and visGP as detailed in the supporting information.

- **visgp_fit.R**: Fits the visGP model.
- **bora.R**: Runs the boraGP method.
- Additional files in the directory are used to compare the resulting model parameters and log-likelihoods from the simulations.
#### cov_comp

This contains **cov_comp.R** which generates the plots for Figure 2 of the paper, using the U-shaped domain.



