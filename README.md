# hmmTMB

This R package implements flexible hidden Markov models, based on Template Model Builder (TMB): flexible state-dependent distributions, transition probability structures, random effects, and smoothing splines. 

## Preprint

The statistical background, as well as details about the implementation of the package, and several example analyses, are presented in the following preprint.

[Michelot, T. (2022). hmmTMB: Hidden Markov models with flexible covariate effects in R. arXiv:2211.14139.](https://arxiv.org/abs/2211.14139)

## Package installation

The package is available on CRAN, and the stable version can therefore be installed using
```
install.packages("hmmTMB")
```

The development version of the package can be installed from Github using devtools,
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

 - ['Analysing time series data with hidden Markov models in hmmTMB'](https://github.com/TheoMichelot/hmmTMB/blob/master/vignettes/hmmTMB_workflow.pdf): Overview of package workflow, using detailed example based on analysis of energy prices. This is a good starting point to learn how to use the package.
 
 - ['Bayesian inference in hmmTMB'](https://github.com/TheoMichelot/hmmTMB/blob/master/vignettes/hmmTMB_example_stan.pdf): Description of workflow for Bayesian analysis in hmmTMB, including specifying priors, and extracting posterior samples.
 
 - ['Advanced features of hmmTMB'](https://github.com/TheoMichelot/hmmTMB/blob/master/vignettes/hmmTMB_advanced_features.pdf): Description of some other useful functionalities, including (semi-)supervised learning, parameter constraints, selection of initial parameter values, etc.
 
 - ['General dependence structures in hmmTMB'](https://github.com/TheoMichelot/hmmTMB/blob/master/vignettes/hmmTMB_general_structures.pdf): Implementation details for hidden Markov models (HMMs) with non-standard dependence structures, including hidden semi-Markov models, higher-order HMMs, autoregressive HMMs, and coupled HMMs.
 
 - ['List of distributions in hmmTMB'](https://github.com/TheoMichelot/hmmTMB/blob/master/vignettes/hmmTMB_dist_list.pdf): List of observation distributions currently available in hmmTMB.
  
 - ['Flexible animal movement modelling using hmmTMB'](https://github.com/TheoMichelot/hmmTMB/blob/master/vignettes/hmmTMB_example_movement.pdf): Description of wild haggis movement analysis, illustrating how non-parametric covariate effects can be included. This includes two different types of movement models: (1) correlated random walks based on step lengths and turning angles, and (2) correlated random walks based on locations directly.
 
 - ['Occupancy modelling using hmmTMB'](https://github.com/TheoMichelot/hmmTMB/blob/master/vignettes/hmmTMB_example_occupancy.pdf): Analysis of occupancy data set of crossbill from [Kéry et al. (2013)](https://onlinelibrary.wiley.com/doi/10.1111/jbi.12087).
