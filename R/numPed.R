################################################
#Adapted from part of the 'inverseA' function
# written by Jarrod Hadfield
#in the 'MCMCglmm' package
################################################

numPed<-function(pedigree)
{
      
  if(length(which(pedigree[, 2] == 0)) > 0) pedigree[which(pedigree[, 2] == 0), 2] <- NA
  if(length(which(pedigree[, 3] == 0)) > 0) pedigree[which(pedigree[, 3] == 0), 3] <- NA


  if (all(is.na(pedigree[, 2])) & all(is.na(pedigree[, 3]))) {
            stop("All dams and sires are missing")
      }
  if (dim(pedigree)[2] != 3) {
            stop("pedigree must have three columns: ID, Dam and Sire")
        }
  if (sum((na.omit(pedigree[, 2]) %in% pedigree[, 1]) == FALSE) > 0 & any(is.na(pedigree[, 2]) == FALSE)) {
            stop("individuals appearing as dams but not in pedigree")
        }
  if (sum((na.omit(pedigree[, 3]) %in% pedigree[, 1]) == FALSE) > 0 & any(is.na(pedigree[, 3]) == FALSE)) {
            stop("individuals appearing as sires but not in pedigree")
        }
  if (sum(duplicated(pedigree[, 1])) > 0) {
            stop("some individuals appear more than once in the pedigree")
        }
  numeric.pedigree <- matrix(as.integer(-998), dim(pedigree)[1], dim(pedigree)[2])
  numeric.pedigree[, 1] <- as.integer(seq(1, dim(pedigree)[1], 1))
  numeric.pedigree[, 2] <- match(pedigree[, 2], pedigree[, 1], nomatch = -998)
  numeric.pedigree[, 3] <- match(pedigree[, 3], pedigree[, 1], nomatch = -998)
  dnmiss <- which(numeric.pedigree[, 2] != -998)
  snmiss <- which(numeric.pedigree[, 3] != -998)
  bnmiss <- which(numeric.pedigree[, 2] != -998 & numeric.pedigree[, 3] != -998)

  if (length(intersect(numeric.pedigree[, 2][dnmiss], numeric.pedigree[, 3][snmiss])) > 0 & (length(dnmiss) > 0) & (length(snmiss) > 0)) {
            warning("Dams appearing as Sires")
        }
  if (any(numeric.pedigree[, 2][dnmiss] > numeric.pedigree[, 1][dnmiss]) & (length(dnmiss) > 0)) {
            stop("Dams appearing before their offspring: first use 'fixPedigree' from pedantics")
        }
  if (any(numeric.pedigree[, 3][snmiss] > numeric.pedigree[, 1][snmiss]) & (length(snmiss) > 0)) {
            stop("sires appearing before their offspring: first use 'fixPedigree' from pedantics")
        }

numeric.pedigree
}


