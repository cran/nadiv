grfx <- function(n, G, incidence = NULL, saveIncidence = FALSE, output = "matrix"){

  d <- dim(G)[1]
  if(all(G == G[1,1]) & d > 1) warning("variance-covariance matrix 'G' may have caused 'chol.default(G)' error.  If so, consider subtracting 0.0001 from the covariances to make correlations < 1 or >-1")
  Mg <- as(chol(G), "dtCMatrix")
  if(is.null(incidence)){
     if(any(ls(envir = globalenv() ) == "nadiv_prev_Mincidence")){
       warning("using previous incidence matrix")
     } else{
          if(saveIncidence){
	     nadiv_prev_Mincidence <<- chol(as(diag(1, n, n), "dgCMatrix"))
             } else{
                  nadiv_prev_Mincidence <- chol(as(diag(1, n, n), "dgCMatrix"))
               }
          warning("Incidence matrix used = Identity matrix")
       }
  } else{
       if(saveIncidence){
          nadiv_prev_Mincidence <<- chol(incidence)
       } else{
            nadiv_prev_Mincidence <- chol(incidence)
         }
    }
 #Must be "dgCMatrix" # chol(incidence), e.g., chol(as(diag(1, n, n), "dgCMatrix"))

  M <- kronecker(nadiv_prev_Mincidence, Mg)
  Z <- Matrix(rnorm(n*d), nrow = 1)
  X <- Matrix((Z %*% M)@x, ncol = d, byrow = TRUE)

 return(as(X, output))
}

