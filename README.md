# hmmTMB

This R package implements flexible hidden Markov models, based on Template Model Builder (TMB): flexible state-dependent distributions, transition probability structures, random effects, and smoothing splines.

## Package installation

Install the package from Github using devtools,
```
devtools::install_github("TheoMichelot/hmmTMB")
```

## Package documentation

To find help files for the methods implemented in the package, search for help using the name of the corresponding class, e.g.,
```
?MarkovChain
?Observation
?HMM
```

We describe functionalities of the package in several vignettes:

 - ['Analysing time series data with hidden Markov models in hmmTMB'](https://github.com/TheoMichelot/hmmTMB/raw/master/vignettes/hmmTMB_workflow.pdf): Overview of package workflow, using detailed example based on analysis of energy prices.
 
 - ['Flexible animal movement modelling using hmmTMB'](https://github.com/TheoMichelot/hmmTMB/raw/master/vignettes/hmmTMB_example_movement.pdf): Description of wild haggis movement analysis, illustrating how non-parametric covariate effects can be included.
 
 - ['Occupancy modelling using hmmTMB'](https://github.com/TheoMichelot/hmmTMB/raw/master/vignettes/hmmTMB_example_occupancy.pdf): Analysis of occupancy data set of crossbill from [Kéry et al. (2013)](https://onlinelibrary.wiley.com/doi/abs/10.1111/jbi.12087).
 
 - ['Bayesian inference in hmmTMB'](https://github.com/TheoMichelot/hmmTMB/raw/master/vignettes/hmmTMB_example_stan.pdf): Description of workflow for Bayesian analysis in hmmTMB, including specifying priors, and extracting posterior samples.
 
 - ['Advanced features of hmmTMB'](https://github.com/TheoMichelot/hmmTMB/raw/master/vignettes/hmmTMB_advanced_features.pdf): Description of some other useful functionalities, including (semi-)supervised learning, parameter constraints, selection of initial parameter values, etc.
 