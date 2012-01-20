aiFun <- function(model = NULL, AI.vec = NULL, inverse = TRUE, Dimnames=NULL)
{
  if(!is.null(model)){
    AI.vec <- model$ai
  }

  AI.cov <- Tri2M(AI.vec, lower.tri = FALSE, reverse = TRUE) 
    if(inverse == FALSE) AI.cov <- solve(AI.cov)    
 
  AI.cor <- cov2cor(AI.cov)
  AI.mat <- matrix(0, nrow=dim(AI.cov)[1], ncol=dim(AI.cov)[2])
  AI.mat[upper.tri(AI.cor, diag=FALSE)] <- AI.cor[upper.tri(AI.cor, diag=FALSE)]
  AI.mat[lower.tri(AI.cov, diag=TRUE)] <- AI.cov[lower.tri(AI.cov, diag=TRUE)]
  if(is.null(Dimnames)){Dimnames <- names(model$gammas)}
  dimnames(AI.mat) <- list(Dimnames, Dimnames)
return(round(AI.mat, 5))
}

