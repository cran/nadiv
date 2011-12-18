makeDsim <- function(pedigree, N, invertD = FALSE, calcSE = FALSE)
{
  numeric.pedigree <- numPed(pedigree)
  n <- dim(pedigree)[1]
  alleles <- matrix(-998, nrow = n, ncol=2) 
  dfounders <- which(numeric.pedigree[, 2] == -998)
  sfounders <- which(numeric.pedigree[, 3] == -998)
  uniqp <- c(unique(numeric.pedigree[, 2])[-1], unique(numeric.pedigree[, 3])[-1])
unique(numeric.pedigree[, 3])[-1]
  ndfounders <- length(dfounders)
  nsfounders <- length(sfounders)

  alleles[dfounders, 1] <- as.integer(seq(1, ndfounders, 1)) 
  alleles[sfounders, 2] <- as.integer(seq(ndfounders+1, (ndfounders + nsfounders), 1))
  nonfounders <-  which(alleles == -998)
  segregation <- matrix(ceiling(runif(length(nonfounders)*N, min = 0, max = 2)),ncol=N)
  dalleles <- matrix(alleles[,1], nrow = dim(alleles)[1], ncol = N)
  salleles <- matrix(alleles[,2], nrow = dim(alleles)[1], ncol = N)

  approxD.tmp <- makeD(pedigree, invertD = invertD)
  approxD <- sm2list(approxD.tmp$D, colnames = c("row", "column", "D"))
  D.nonself <- which(approxD[,3] < 1) 
  approxD.nonself <- approxD[D.nonself,]
  simD.nonself <- matrix(NA, nrow = dim(approxD.nonself)[1], ncol = N) 

print("starting to make Dsim")
    start <- 1
    order.ind <- order(uniqp)
    uniqp <- uniqp[order.ind]
    for(j in uniqp){
      mothered <- NULL
      fathered <- NULL
      if(any(numeric.pedigree[,2] == j)) mothered <- which(numeric.pedigree[,2] == j)
      if(any(numeric.pedigree[,3] == j)) fathered <- which(numeric.pedigree[,3] == j)
      jalleles <- rbind(dalleles[j, ], salleles[j, ])
      lmom <- length(mothered)
      lfath <- length(fathered)
      segregation.tmp <- matrix(segregation[start:(start - 1 + lmom + lfath), ], ncol = N)
      jgametes <- matrix(vapply(1:(lmom+lfath), FUN = vchoose, FUN.VALUE = vector("numeric", N), N = N, which2 = segregation.tmp, values2 = jalleles), ncol = N, byrow = TRUE) 
      start <- start + lmom + lfath 
     if(!is.null(mothered)) dalleles[mothered, ] <- jgametes[1:lmom, ]
     if(!is.null(fathered)) salleles[fathered, ] <- jgametes[(lmom + 1):lfath, ]
    }

  stack.alleles<-rbind(dalleles, salleles)
  simD.nonself <- apply(stack.alleles, MAR = 2, FUN = IBD2, pairs = approxD.nonself[,1:2], n=n) 
  approxD.nonself$simD <- apply(simD.nonself, MAR = 1, FUN = mean)
  listDsim <- NULL
  if(calcSE) {
     approxD.nonself$Dse <- vapply(approxD.nonself$simD, FUN = SEfun, FUN.VALUE = vector("numeric", 1), N)
     listDsim <- approxD.nonself
   } 

  Dsim.row<-c(approxD.nonself[,1], 1:n)
  Dsim.col<-c(approxD.nonself[,2], 1:n)
  Dsim.x<-c(approxD.nonself[,4], rep(1, n))
  order.index<-order(Dsim.col + Dsim.row/(n+1), decreasing=FALSE)
  Dsim<-Matrix(0, n, n)
  Dsim@uplo<-"L"
  Dsim@i<-as.integer(Dsim.row[order.index]-1)
  Dsim@p<-as.integer(c(match(1:n, Dsim.col[order.index]), length(order.index)+1)-1)
  Dsim@x<-Dsim.x[order.index]
print("Dsim made")
  logDetDsim <- determinant(Dsim, logarithm = TRUE)$modulus[1]

 
  if(invertD){
    print("starting to invert Dsim")
    Dsiminv<-solve(Dsim)
    Dsiminv@Dimnames <- list(pedigree[,1], NULL)
    print("done inverting Dsim")
    listDsiminv<-sm2list(Dsiminv, rownames=pedigree[,1], colnames=c("row", "column", "simDinverse"))
    Dsim <- as(Dsim, "dgCMatrix")
 return(list(A=approxD.tmp$A, D=approxD.tmp$D, logDetD = approxD.tmp$logDet, Dinv=approxD.tmp$Dinv, listDinv=approxD.tmp$listDinv, Dsim=Dsim, logDetDsim = logDetDsim, Dsiminv=Dsiminv, listDsim=listDsim, listDsiminv=listDsiminv))
  } else{
    return(list(A=approxD.tmp$A, D=approxD.tmp$D, logDetD = approxD.tmp$logDet, Dsim=Dsim, logDetDsim = logDetDsim, listDsim=listDsim))
    } 

}


