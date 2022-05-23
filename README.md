# NCICS: Numerical considerations for ICS

This repository contains the code to reproduce the artificial examples of the 
following paper: Numerical considerations and a new implementation for ICS, 
Archimbaud A., Drmac Z., Nordhausen K., Radojicic U. and A. Ruiz-Gazen, 2022 (submitted).


## Organization of the repository

### Data folder

Two new anonymized data sets, HTP2 and HTP3 are available in the *data* folder in `.rda` format. They are described in the paper and in the *R/data.R* file.


### Figures folder

This folder is the supplementary material of the paper, containing all the plots
included in the paper as well as additional ones.

### R folder

It contains three R scripts:

- `all_applications.R`: the script to reproduce all the results and figures from the 
paper.

- `data.R`: the description of the new data sets, HTP2 and HTP3.

- `ICS_utilities.R`: all self-written R functions to perform the analysis. The
main two are: `ICS_QR_full_rank()` which implements the ICSQR algorithm described in
the paper and `ICS_QR_not_full_rank()` in case the data are not full rank and we want
to perform first a dimension reduction based on a pivoted QR decomposition.