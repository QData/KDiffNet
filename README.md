# KDiffNet

This repository contains code for estimating Differential Networks using extra evidence as described in the following paper:  
>[Beyond Data Samples: Aligning Differential Networks Estimation with Scientific Knowledge](https://arxiv.org/abs/2004.11494)  
> Arshdeep Sekhon, Zhe Wang, Yanjun Qi
> International Conference on Artificial Intelligence and Statistics (AISTATS), 2022



## Estimating Differential Matrix using KDiffNet 

To estimate the differential matrix `estimated_delta` given two sets of samples, `X_c` and `X_d` and edge level known knowledge `W_E` and node group knowledge `G_V`:
```{r}
sourceDirectory("kdiffnet/")
estimated_delta = kdiffnet(X_c, X_d, W_E, G_V, epsilon, lambda)
```
Here, `lambda` and `epsilon` are hyperparameters. 

## Generating Simulated Data 

To generate a ground truth differential matrix `truth_delta`, we first need to specify its generating parameters. To generate `truth_delta` from the edge level matrix `W_E` and known groups `G_V`, with a sparsity level of `delta_sparsity` and `individual_sparsity`. 
```{r}
graphs = simulate_w_group(W, g, delta_sparsity, individual_sparsity)
Omega_d = graphs[[2]]
Omega_c = graphs[[1]]
truth_delta=(Omega_d - Omega_c)
```

We can then generate two sets of samples, `X_c` and `X_d`, each following a Gaussian Distribution with inverse covariance matrix as Omega_c` and `Omega_d ` respectively:

```{r}
f <- function(m) class(try(solve(m),silent=T))=="matrix"
library(Hmisc)
library(R.utils)
library(boot)
library(MASS)
covcTrue =  solve(Omega_c)
covdTrue =  solve(Omega_d)
X_c = mvrnorm(nc, mu = rep(0,nrow(W_E)), covcTrue)
X_d = mvrnorm(nd, mu = rep(0,nrow(W_E)), covdTrue)
```
Here, `nc` and `nd` are the number of samples we wish to generate. 
