aiFun <- function(model = NULL, AI.vec = NULL, inverse = TRUE, Dimnames=NULL)
{
  if(!is.null(model)){
    AI.vec <- model$ai
  }

  dimAI <- sqrt(length(AI.vec) * 2 + 0.25) - 0.5
  AI.cov <- matrix(0, dimAI, dimAI)
  AI.cov[which(upper.tri(AI.cov, diag = TRUE) == TRUE)] <- AI.vec
  AI.cov[which(lower.tri(AI.cov) == TRUE)]<-t(AI.cov)[which(lower.tri(AI.cov) == TRUE)]
    if(inverse == FALSE) AI.cov <- solve(AI.cov)    
 
  AI.cor <- cov2cor(AI.cov)
  AI.mat <- matrix(0, nrow=dim(AI.cov)[1], ncol=dim(AI.cov)[2])
  AI.mat[upper.tri(AI.cor, diag=FALSE)] <- AI.cor[upper.tri(AI.cor, diag=FALSE)]
  AI.mat[lower.tri(AI.cov, diag=TRUE)] <- AI.cov[lower.tri(AI.cov, diag=TRUE)]
  if(is.null(Dimnames)){Dimnames <- names(model$gammas)}
  dimnames(AI.mat) <- list(Dimnames, Dimnames)
       

return(AI = round(AI.mat, 5))
}

