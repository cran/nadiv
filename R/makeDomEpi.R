makeDomEpi <- function(pedigree, output = c("AD", "DD", "both"), Dinverse = FALSE)
{
  type <- match.arg(output)
  Dout <- makeD(pedigree, invertD=Dinverse)

  if(type == "AD"){
    AD <- Dout$A * Dout$D
    logDetAD <- determinant(AD, logarithm = TRUE)$modulus[1]
    ADinv <- Matrix(solve(AD), sparse=TRUE)
    listADinv <-sm2list(ADinv, rownames=pedigree[,1], colnames=c("row", "column", "ADinverse"))
    ADinv <- as(ADinv, "dgCMatrix")
    DD <- NULL
    logDetDD <- NULL
    DDinv <- NULL
    listDDinv <- NULL
    }      

  if(type == "DD"){
    DD <- Dout$D * Dout$D
    logDetDD <- determinant(DD, logarithm = TRUE)$modulus[1]
    DDinv.tmp <- Matrix(solve(DD), sparse=TRUE)
    listDDinv<-sm2list(DDinv.tmp, rownames=pedigree[,1], colnames=c("row", "column", "DDinverse"))
    DDinv <- as(DDinv.tmp, "dgCMatrix")
    AD <- NULL
    logDetAD <- NULL
    ADinv <- NULL
    listADinv <- NULL
   }

  if(type == "both"){
    AD <- Dout$A * Dout$D
    logDetAD <- determinant(AD, logarithm = TRUE)$modulus[1]
    ADinv <- Matrix(solve(AD), sparse=TRUE)
    listADinv <-sm2list(ADinv, rownames=pedigree[,1], colnames=c("row", "column", "ADinverse"))
    ADinv <- as(ADinv, "dgCMatrix")
    DD <- Dout$D * Dout$D
    logDetDD <- determinant(DD, logarithm = TRUE)$modulus[1]
    DDinv.tmp <- Matrix(solve(DD), sparse=TRUE)
    listDDinv<-sm2list(DDinv.tmp, rownames=pedigree[,1], colnames=c("row", "column", "DDinverse"))
    DDinv <- as(DDinv.tmp, "dgCMatrix")
    }
return(list(D=Dout$D, logDetD = Dout$logDet, AD=AD, logDetAD = logDetAD, DD=DD, logDetDD = logDetDD, Dinv=Dout$Dinv, ADinv=ADinv, DDinv=DDinv, listDinv=Dout$listDinv, listADinv=listADinv, listDDinv=listDDinv))

}

