IBD2 <- function(pairs, genos, n){
  apply(pairs, FUN = IBD, MARGIN = 1, genotypes = genos, n = n)
}

