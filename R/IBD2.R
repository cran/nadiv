IBD2 <- function(pairs, genos, n){
  apply(pairs, FUN = IBD, MAR = 1, genotypes = genos, n = n)
}

