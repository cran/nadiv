drfx <- function(G, fac, dataf){
    dataf[, fac] <- as.factor(dataf[, fac])
    d <- dim(G)[1]
    if(all(G == G[1,1]) & d > 1) warning("variance-covariance matrix 'G' may have caused 'chol.default(G)' error.  If so, consider subtracting 0.0001 from the covariances to make correlations < 1 or >-1")

   Z <- sparse.model.matrix(as.formula(paste("~", fac, " - 1", sep = "")), dataf)
   M <- suppressWarnings(grfx(n = dim(Z)[2], G = G))
   fx <- matrix(NA, nrow = dim(dataf)[1], ncol = d)
   fx[as.integer(dimnames(Z)[[1]]), ] <- sapply(seq.int(d), FUN = function(c){ (Z %*% M[, c])@x})


return(list(fx = fx, Z = Z))
}

