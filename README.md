# BayesianFunctionalAutoregression
This repository contains the code and materials developed for my thesis on Bayesian modelling for autoregressive functional data
## Introduction
Presentation of the problem and actual state-of-the-art solutions. Discussion on the advantages of adopting a Bayesian framework in this context. 

## Statistical Framework
Overview of the main statistical tools used throughout the work: 
- Hierarchical modelling and posterior conjugacy. 
- Gibbs sampling, Metropolis Hastings, and Multiple Try MH algorithms.
- Empirical Bayes for prior inizializations. 
- Spike and Slab prior for variable selection.
- Theoretical background for Fourier and splines decompositions. 

## Bayesian Modelling
Development fo the Bayesian functional autoregressive model, including: 
- Derivation of conjugacy and construction of the Gibbs Sampler.
- Extension to models with Phi decomposition.
- Modelling of no autoregressive, first-order autoregressive, and higher-order autoregressive effects with variable selection. 


## Synthetic Study
Simulation of random autoregressive curves under different scenarios and performance measuring of the model in component identification and parameter recovery. 

## Parkinson Data application
Application of the model proposed to real world data: functional measurements of hand tremors in Parkinson’s disease patients
