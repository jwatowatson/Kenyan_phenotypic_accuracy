**Improving statistical power in severe malaria genetic association studies by augmenting phenotypic precision** (eLife in press, 2021)

The latest preprint version can be found at: https://www.biorxiv.org/content/10.1101/2021.04.16.440107v2


# Overview of repository

The underlying work is broken down into a few RMarkdown and R scripts:

* Data_Prep: this puts together the analysis datasets into a manageable set. This script can only be run with access to the original data files - this needs permission via data access committees (see details below)
* Mixture_modelling: this is the main script for the diagnostic model of severe malaria. This generates Figures 1 and 3 of the paper. This script can be run using the curated dataset provided. Note that the Bayesian fitting for the diagnostic model takes a long time on a standard computer (approx 10 hours)
* Direct_typed_Association_Study: association study using the directly typed SNPs. This can only be run with access to the directly typed polymorphism data
* Extra_analyses: a few extra bits using the estimated probabilities of severe malaria, using data that are not open access (positive blood cultures, haematocrits)
* Simulation_study_weightedLikelihood: this implements the simulations given in Appendix 12 of the paper, showing how the weighted likelihood works


Any questions contact me: jwatowatson at gmail dot com


## Software packages needed

To run the models you will need a few widely used R packages plus two key model fitting packages:

* *rstan* (run the Bayesian models)
* *mgcv* (fit the generalised additive models using penalised splines)
