makeD <- function(pedigree, invertD=TRUE){
  numeric.pedigree <- numPed(pedigree) 
  numeric.pedigree[numeric.pedigree == -998] <- NA
  A <- makeA(numeric.pedigree)
  listA <- sm2list(A, rownames=numeric.pedigree[,1], colnames=c("row", "column", "A"))

  noself.listA <- listA[!listA[, 1] ==listA[, 2], ]
  exclude <- which(is.na(numeric.pedigree[, 2]) | is.na(numeric.pedigree[, 3]))
  tmp.listA <- noself.listA[!noself.listA[,1] %in% exclude & !noself.listA[,2] %in% exclude,]

  iparents <- matrix(c(numeric.pedigree[tmp.listA[,1], 2],numeric.pedigree[tmp.listA[,1], 3],numeric.pedigree[tmp.listA[,2], 2],numeric.pedigree[tmp.listA[,2], 3]), nrow=4, ncol=dim(tmp.listA)[1], byrow=TRUE)

  increment <- round(((ceiling(dim(iparents)[2]*0.05)-(dim(iparents)[2]*0.05))*20)+0.1)
  if(!increment==0) {
    addition <- matrix(1, nrow=4, ncol=increment)
    iparents <- cbind(iparents, addition)
  }

  iparents_r <- dim(iparents)[2]
  N <- dim(A)[1]
  dim_aij <- dim(listA)[1]
  rcA <- t(listA[, 1:2])
  aij <- listA[, 3]/2
  answer <- as.single(rep(0,iparents_r))
  cat(paste("starting to make D..."))


  For_out <- .Fortran(dij, iparents=as.integer(iparents), iparents_r=as.integer(iparents_r), dim_aij=as.integer(dim_aij), rcA=as.integer(rcA), aij=as.single(aij), N=as.integer(N), answer=as.single(answer))$answer


  if(!increment==0){
    Dijs <- For_out[-c((length(For_out)-increment+1):length(For_out))]
    } else{Dijs <- For_out}
  cat(paste(".done", "\n"))
	
  tmp.listD <- data.frame(Row=tmp.listA[,1], Column=tmp.listA[,2], D=Dijs)
  tmp.listD2 <- tmp.listD[which(!tmp.listD[,3]==0), ]

  D.row <- c(tmp.listD2[,1], 1:N)
  D.col <- c(tmp.listD2[,2], 1:N)
  D.x <- c(tmp.listD2[,3], rep(1, N))
  order.index <- order(D.col + D.row/(N+1), decreasing=FALSE)
  D <- Matrix(0, N, N)
  D@uplo <- "L"
  D@i <- as.integer(D.row[order.index]-1)
  D@p <- as.integer(c(match(1:N, D.col[order.index]), length(order.index)+1)-1)
  D@x <- D.x[order.index]
  
  logDet <- determinant(D, logarithm = TRUE)$modulus[1]
 
  if(invertD){
    Dinv <- solve(D)
    Dinv@Dimnames <- list(pedigree[,1], NULL)
    listDinv <- sm2list(Dinv, rownames=pedigree[,1], colnames=c("row", "column", "Dinverse"))
    D <- as(D, "dgCMatrix")
 return(list(A=A, D=D, logDet = logDet, Dinv=Dinv, listDinv=listDinv))
  } else{
    D <- as(D, "dgCMatrix")
    return(list(A=A, D=D, logDet = logDet))
    } 
}


