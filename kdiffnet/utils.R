# additional functions for kdiffnet
compute_cov <- function(X,covType = "cov"){

  result = X
  if (is.data.frame(X)){
    result = as.matrix(X)
  }
  if (!isSymmetric(X)){
    try(if (covType %in% c("cov","kendall") == FALSE) stop("The cov/cor type you specifies is not include in this package. Please use your own function to obtain the list of cov/cor and use them as the input of simule()"))
    if (covType == "cov"){
        result= stats::cov(result)
    }
    if (covType == "kendall"){
      result= cor.fk(result)
    }
  }
  return(result)
}


checkInv <- function(m) class(try(solve(m),silent=T))=="matrix"

softThre <- function(x, lambda){
    result = sign(x) * pmax(abs(x)-lambda, 0)
    return(result)
}

backwardMap <-function(covMatrix, thre = "soft"){
  niuList = 0.001 * (0:1000) * max(covMatrix)
  bestDet = det(softThre(covMatrix, 0.001))
  bestniu = 0.001

  if (thre == "soft"){
    for (i in 1:1000){
      if (bestDet < det(softThre(covMatrix, niuList[i]))){
        bestDet = det(softThre(covMatrix, niuList[i]))
        bestniu = niuList[i]
      }
    }
    return(solve(softThre(covMatrix, bestniu)))
  }

  if (thre == "hard"){
    for (i in 1:1000){
      if (checkInv(hardThre(covMatrix, niuList[i]))){
        bestniu = niuList[i]
        break
      }
    }
    return(solve(hardThre(covMatrix, bestniu)))
  }
}


