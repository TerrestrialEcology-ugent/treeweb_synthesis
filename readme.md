# TREEWEB synthesis 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1462802.svg)](https://doi.org/10.5281/zenodo.1462802)

This repository contain all scripts and data needed to reproduce the analysis and main figures presented in the published paper: "Forest fragmentation modulates effects of tree species richness and composition on ecosystem multifunctionality", full text available [here](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2653) pdf are available upon request. 

## How to

* the R-scripts to reproduce the analysis are in the **scripts** folder, the main file is _synthesis\_script.r_
* code to run the 10-fold cross-validation is in the _kfold\_script.r_, beware running this takes some time (!)
* the Stan code of the models is in the **model** folder
* the shiny app to explore the sensitivity of the effect of tree richness and fragmentation intensity on multifunctional desirability is in the **shiny** folder, in an R session loading the shiny package and running "runApp("shiny/")" should do the charm

## Note on R and packages version

All these analysis were ran with the following R and packages version:

* R v3.4.4
* rstan v2.17
* dplyr v0.7.6


