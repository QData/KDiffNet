# generate simulated graphs
simulate_w_group <- function(W,g,delta_sparsity,individual_delta_sparsity){
library(boot)
W=(W+t(W))/2
p=nrow(W)
delta_w = matrix(0,p,p)
delta_g = matrix(0,p,p)
  for (j in 1:p) {
    for (k in 1:p) {
      if (g[j] == g[k] & g[j]!=0.0){
        delta_g[j,k] = 1/3
        delta_g[k,j] = delta_g[j,k]
      }else{
        delta_g[j,k] = 0.0
        delta_g[k,j] = delta_g[j,k]
      }
    }
  }

for(i in 1:p){
for(j in i:p){
delta_w[i,j]=(1/3)*sum(inv.logit(-W[i,j])>delta_sparsity)
delta_w[j,i]=delta_w[i,j]
}
}
individual = matrix(0,p,p)
for(i in 1:p){
for(j in i:p){
individual[i,j]=(1/3)*rbinom(1,1,individual_delta_sparsity)
individual[j,i]=individual[i,j]
}
}
graphs =  list()
graphs[[1]] = individual

graphs[[2]] = delta_w + delta_g + individual
graphs[[3]] = delta_w + delta_g


for(i in 1:2){
for(j in 1:p){
graphs[[i]][j,j]=1.0
}
}
I = diag(1,p,p)
for (i in 1:2) {
    graphs[[i]] = (graphs[[i]] - I) + (abs(min(eigen(graphs[[i]] - I)$value)) + 1) * I
  }
return(graphs)


}
