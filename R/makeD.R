makeD <- function(pedigree, invertD=TRUE)
{
  Aout <- makeA(pedigree)
  Adense <- Aout$Adense
  A <- Aout$Asparse
  listA <- sm2list(A, rownames=pedigree[,1], colnames=c("row", "column", "A"))

  numeric.pedigree <- numPed(pedigree)
  noself.listA <- listA[!listA[,1]==listA[,2],]
  exclude <- which(numeric.pedigree[,2]==-998 | numeric.pedigree[,3]==-998)
  tmp.listA <- noself.listA[!noself.listA[,1] %in% exclude & !noself.listA[,2] %in% exclude,]

  iparents <- matrix(c(numeric.pedigree[tmp.listA[,1], 2],numeric.pedigree[tmp.listA[,1], 3],numeric.pedigree[tmp.listA[,2], 2],numeric.pedigree[tmp.listA[,2], 3]), nrow=4, ncol=dim(tmp.listA)[1], byrow=TRUE)

  increment <- round(((ceiling(dim(iparents)[2]*0.05)-(dim(iparents)[2]*0.05))*20)+0.1)
  if(!increment==0) {
    addition <- matrix(1, nrow=4, ncol=increment)
    iparents <- cbind(iparents, addition)
  }
  iparents_r <- dim(iparents)[2]
  dim_a <- dim(Adense)[1]
  answer <- as.single(rep(0,iparents_r))
  print("starting to make D")
  For_out <- .Fortran(dij, iparents=as.integer(iparents), iparents_r=as.integer(iparents_r), Adense=as.single(Adense), dim_a=as.integer(dim_a), answer=as.single(answer))$answer

  if(!increment==0){
    Dijs <- For_out[-c((length(For_out)-increment+1):length(For_out))]
    } else{Dijs <- For_out}
  print("D made")
	
  n_ped <- dim(numeric.pedigree)[1]
  tmp.listD <- data.frame(Row=tmp.listA[,1], Column=tmp.listA[,2], D=Dijs)
  tmp.listD2 <- tmp.listD[which(!tmp.listD[,3]==0), ]

  D.row <- c(tmp.listD2[,1], 1:n_ped)
  D.col <- c(tmp.listD2[,2], 1:n_ped)
  D.x <- c(tmp.listD2[,3], rep(1, n_ped))
  order.index <- order(D.col + D.row/(n_ped+1), decreasing=FALSE)
  D <- Matrix(0, n_ped, n_ped)
  D@uplo <- "L"
  D@i <- as.integer(D.row[order.index]-1)
  D@p <- as.integer(c(match(1:n_ped, D.col[order.index]), length(order.index)+1)-1)
  D@x <- D.x[order.index]
  
  logDet <- determinant(D, logarithm = TRUE)$modulus[1]
  A <- as(A, "dgCMatrix")
 
  if(invertD){
    print("starting to invert D")
    Dinv <- solve(D)
    Dinv@Dimnames <- list(pedigree[,1], NULL)
    print("done inverting D")
    listDinv <- sm2list(Dinv, rownames=pedigree[,1], colnames=c("row", "column", "Dinverse"))
    D <- as(D, "dgCMatrix")
 return(list(A=A, D=D, logDet = logDet, Dinv=Dinv, listDinv=listDinv))
  } else{
    D <- as(D, "dgCMatrix")
    return(list(A=A, D=D, logDet = logDet))
    } 
}

