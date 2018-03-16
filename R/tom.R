#' Permutation
perm <- function (n, r, v = 1:n)
{
  if (r == 1) 
    matrix(v, n, 1)
  else if (n == 1) 
    matrix(v, 1, r)
  else {
    X <- NULL
    for (i in 1:n) X <- rbind(X, cbind(v[i], perm(n - 1, r - 1, v[-i])))
    X
  }
}

#' Topology overlap matrix
TOMsimilarity <- function(adj){
  diag(adj) <- rep(0, nrow(adj))
  L <- adj %*% adj
  kij <- colSums(adj)
  K <- array(0L, dim(adj))
  
  p <- perm(nrow(adj), 2)
  for (idx in 1:nrow(p)){
    K[p[idx,1], p[idx,2]] <- min(kij[p[idx,1]], kij[p[idx, 2]])
  }
  W <- (L + adj)/(K + 1 - adj)
  diag(W) <- rep(1, nrow(W))
  
  return(W)
  
}

#' Adjacenty matrix
#' 
#' @importFrom stats cor
adjacency <- function(datExpr, type = "signed", power = 6){
  sim <- cor(datExpr)
  if (type == "signed"){
    return(abs(0.5 + 0.5*sim)^power)
  }
  else if (type == "unsigned"){
    return(abs(sim)^power)
  }
  else{
    stop("type should be one of 'signed' or 'unsigned'")
  }
  
}