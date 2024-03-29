# NCICS: Numerical considerations for ICS

This repository contains the code to reproduce the examples of the 
following paper: Numerical considerations and a new implementation for ICS, 
Archimbaud, A., Drmac, Z., Nordhausen, K., Radojicic, U. & Ruiz-Gazen, A. (2023). SIAM Journal on Mathematics of Data Science (SIMODS), Vol.5(1):97–121.
https://epubs.siam.org/doi/10.1137/22M1498759


## Organization of the repository

### Data folder

Two new anonymized data sets, HTP2 and HTP3 are available in the *data* folder in `.rda` format. 
They are described in the paper and in the *R/data.R* file.


### Figures folder

This folder is the supplementary material of the paper, containing all the plots
included in the paper as well as additional ones.

### R folder

It contains three R scripts:

- `all_applications.R`: the script to reproduce all the results and figures from the 
paper.

- `data.R`: the description of the new data sets, HTP2 and HTP3.

- `ICS_utilities.R`: all self-written R functions to perform the analysis. The
 two main functions are: `ICS_QR_full_rank()` which implements the ICSQR algorithm described in
the paper and `ICS_QR_not_full_rank()` in case the data are not full rank and we want
to perform first a dimension reduction based on a pivoted QR decomposition.

### renv.lock

A file *renv.lock* is provided for reproducible purposes. 
