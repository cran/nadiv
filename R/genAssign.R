genAssign <- function(pedigree)
{ 
   n <- dim(pedigree)[1]
   generation <- vector("integer", length = n)
   if(!all(apply(pedigree[, 1:3], MARGIN = 2, FUN = is.numeric)) | any(apply(pedigree[, 1:3], MARGIN = 2, FUN = is.na))){
      pedigree[, 1:3] <- numPed(pedigree[, 1:3])
   }

   Cout <- .C("ga",
	as.integer(pedigree[, 2] - 1),
	as.integer(pedigree[, 3] - 1),
        generation,
	as.integer(n))

  Cout[[3]]
}
