# this script can be run using "Rscript example_kdiffnet.R"
# This script shows one example setting of KDiffNet, the user can modify the different simulation parameters as well as 
# the different estimation parameters
# Step 1: Simulates the ground truth
# Step 2: Estimates using KDiffNet and compares to the ground truth 


setwd("./")
f <- function(m) class(try(solve(m),silent=T))=="matrix"
library(Hmisc)
library(R.utils)
library(boot)
library(MASS)
set.seed(37)
sourceDirectory("kdiffnet/")




#using abide atlas 
brain_dist<-("kdiffnet/abide_distance_matrix.csv")
W=read.csv(brain_dist,sep=',',header=FALSE)
W=as.matrix((W + t(W)) /2)

p=ncol(W)

# """""""""""""""""""""""""STEP 1: Simulate Data Samples""""""""""""""""""""""""""""""
# simulation parameters
# group size
size=10
# proportion of samples to generate: nsamples = rc*p // rd *p
rc=1
rd=1
individual_sparsity = 0.5
# level of sparsity of underlying true delta 
d_index=5
known_index=1
ngroups_list=seq(1,floor(p/size))
ngroups=floor(median(ngroups_list))
# threshold different levels of known W
delta_sparsity_values=quantile(inv.logit(-W), c(.125, .25, .375,.5,.625,.75,.875))
delta_sparsity=delta_sparsity_values[d_index]
WKnown=W
WKnown[which(inv.logit(-W)>delta_sparsity)]=1.0
g=rep(0,p)
# number of nodes in some group
gs = ngroups*size
gid=1
for(g_index in 1:(gs)){
g[g_index]=gid
if(g_index%%10==0){
gid=gid+1
}
}
graphs = simulate_w_group(W,g,delta_sparsity,individual_sparsity)
# Make difference graph
truth_delta=(graphs[[2]]-graphs[[1]])
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""





# Generate samples 


nc=floor(p * rc)
nd=floor(p * rc)
covcTrue =  solve(graphs[[1]])
covdTrue =  solve(graphs[[2]])



# generate training data 
datac = mvrnorm(nc, mu = rep(0,nrow(W)), covcTrue)
datad = mvrnorm(nd, mu = rep(0,nrow(W)), covdTrue)


# generate a separate set of validation data  
valid_datac = mvrnorm(nc, mu = rep(0,nrow(W)), covcTrue)
valid_datad = mvrnorm(nd, mu = rep(0,nrow(W)), covdTrue) 



# SINGLE USE CASE""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#epsilon=1.0
#lambda = 0.01*1.0*sqrt(log(p)/min(nc,nd))
#estimated_delta = kdiffnet(datac, datad, WKnown, g, epsilon, lambda)
#f1_score = evaluate(truth_delta,as.matrix(estimated_delta))
# print(f1_score)






# """""""""""""""""""""""""""""""""""""""""STEP 2: estimate difference by validating over hyperparameter lambda """"""""""""""""""""""""""""""""""""""""""""
# in this example, we only iterate over 10 lambda values
epsilon=1.0
lambda_vals=seq(1:10)
best_f1_score=-1.0
# ------can precompute backmap to speedup computation------------
#Back = precompute_backmap(datac,datad)

for (lambda_c in lambda_vals){
lambda = 0.01*lambda_c*sqrt(log(p)/min(nc,nd))

#estimated_delta = kdiffnet(B=Back, W=WKnown, g=g, epsilon=epsilon, lambda=lambda,precompute=TRUE)
#f1_score = evaluate(truth_delta,as.matrix(estimated_delta))
estimated_delta = kdiffnet(datac, datad, WKnown, g, epsilon, lambda)
f1_score = evaluate(truth_delta,as.matrix(estimated_delta))
#print(paste0(' estiamted training f1 for difference graph: ', f1_score))

if(f1_score>best_f1_score){
best_lambda = lambda
best_f1_score = f1_score
}
}



estimated_delta = kdiffnet(valid_datac, valid_datad, WKnown, g, epsilon, best_lambda)
f1_score = evaluate(truth_delta,as.matrix(estimated_delta))
print(paste0(' f1 score for  estimated difference graph vs truth simulation: ', f1_score))






