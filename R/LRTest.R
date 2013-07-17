LRTest <- function(full, reduced, df = 1){
    lambda <- (-2*(full - reduced))
    lrtP <- pchisq(lambda, df=df, lower.tail = FALSE)
 return(list(lambda = lambda, Pval = lrtP))
}

