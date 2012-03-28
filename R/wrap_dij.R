wrap_dij <- function(x, parents, relatedness, N){
    tmp_parents <- parents[, min(x):max(x)]
    parents_r <- dim(tmp_parents)[2]
    increment <- round(((ceiling(parents_r*0.05)-(parents_r*0.05))*20)+0.1)
    if(!increment == 0){
      addition <- matrix(1, nrow = 4, ncol = increment)
      tmp_parents <- cbind(tmp_parents, addition)
      parents_r <- parents_r + increment
    }

    dim_aij <- dim(relatedness)[1]
    rcA <- t(relatedness[, 1:2])
    aij <- relatedness[, 3]/2
    answer <- as.single(rep(0,parents_r))

    For_out <- .Fortran(dij, iparents=as.integer(tmp_parents), parents_r=as.integer(parents_r), dim_aij=as.integer(dim_aij), rcA=as.integer(rcA), aij=as.single(aij), N=as.integer(N), answer=as.single(answer))$answer
 
    if(!increment == 0){
       dij_out <- For_out[-c((length(For_out) - increment + 1):length(For_out))]
    } else dij_out <- For_out
  dij_out
}

