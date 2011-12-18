makeAA <- function(pedigree)
{
  A <- makeA(pedigree)$Asparse
  AA <- A*A
  logDet <- determinant(AA, logarithm = TRUE)$modulus[1]
  AAinv <- as(solve(AA), "dgCMatrix")
  listAAinv <- sm2list(AAinv, rownames=pedigree[,1], colnames=c("row", "column", "AAinverse"))
  AA <- as(AA, "dgCMatrix")
return(list(AA=AA, logDet = logDet, AAinv=AAinv, listAAinv=listAAinv))    
}

