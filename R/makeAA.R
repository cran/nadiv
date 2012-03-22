makeAA <- function(pedigree)
{
  if(dim(pedigree)[1] < 3000) meth <- "fast" else meth <- "flowmem"
  A <- makeA(pedigree, method = meth)
  AA <- A*A
  logDet <- determinant(AA, logarithm = TRUE)$modulus[1]
  AAinv <- as(solve(AA), "dgCMatrix")
  listAAinv <- sm2list(AAinv, rownames=pedigree[,1], colnames=c("row", "column", "AAinverse"))
  AA <- as(AA, "dgCMatrix")
return(list(AA=AA, logDet = logDet, AAinv=AAinv, listAAinv=listAAinv))    
}

