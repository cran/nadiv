LRTest<-function(full, reduced, df = 1){
    lambda<-(-2*(full - reduced))
    lrtP<-pchisq(lambda, df=df, lower.tail = FALSE, log.p = TRUE)
 list(lambda = lambda, Pval = lrtP)
}

