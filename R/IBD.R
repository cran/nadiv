IBD <- function(pair, genotypes, n){
  rgeno <- c(genotypes[pair[[1]]], genotypes[pair[[1]]+n])
  rgeno.alt <- c(genotypes[pair[[1]]+n], genotypes[pair[[1]]])
  cgeno <- c(genotypes[pair[[2]]], genotypes[pair[[2]]+n])

  if(identical(rgeno, cgeno) | identical(rgeno.alt, cgeno)) return(1)
    else return(0)
}

