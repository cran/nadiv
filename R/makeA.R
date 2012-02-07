makeA <- function(pedigree)
{

  numeric.pedigree <- numPed(pedigree)
  N <- dim(numeric.pedigree)[1]

  tmp.ped <- data.frame(sire = numeric.pedigree[,3], dam = numeric.pedigree[,2], label = as.character(numeric.pedigree[,1]))
  tmp.ped[tmp.ped == -998] <- NA
  ped <- pedigree(tmp.ped[,1], tmp.ped[,2], tmp.ped[,3])
 
  tL <- relfactor(ped, ped@label)
  A <- as(crossprod(tL), "dgCMatrix")
  
A
}

