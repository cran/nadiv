makeD <- function(pedigree, parallel = FALSE, ncores = getOption("cores"), invertD = TRUE, Amethod = NULL, returnA = FALSE)
{

  numeric.pedigree <- numPed(pedigree) 
  numeric.pedigree[numeric.pedigree == -998] <- NA
  N <- dim(pedigree)[1]
  if(is.null(Amethod)){
    if(N < 3000) meth <- "fast" else meth <- "flowmem"
  } else meth <- Amethod
  A <- makeA(numeric.pedigree, method = meth)
  listA <- sm2list(A, rownames=numeric.pedigree[,1], colnames=c("row", "column", "A"))
     if(!returnA) A <- NULL
  listA <- listA[which(listA[,3] != 0), ]
  noself.listA <- listA[!listA[, 1] ==listA[, 2], ]
  exclude <- which(is.na(numeric.pedigree[, 2]) | is.na(numeric.pedigree[, 3]))
  tmp.listA <- noself.listA[!noself.listA[,1] %in% exclude & !noself.listA[,2] %in% exclude,]

  iparents <- matrix(c(numeric.pedigree[tmp.listA[,1], 2],numeric.pedigree[tmp.listA[,1], 3],numeric.pedigree[tmp.listA[,2], 2],numeric.pedigree[tmp.listA[,2], 3]), nrow=4, ncol=dim(tmp.listA)[1], byrow=TRUE)

  if(parallel & dim(iparents)[2] < 80){
     cat(paste("Warning: pedigree too small - 'parallel' set to FALSE instead", "\n"))
     parallel <- FALSE
  }

  if(parallel == FALSE){
     increment <- round(((ceiling(dim(iparents)[2]*0.05)-(dim(iparents)[2]*0.05))*20)+0.1)
     if(!increment==0) {
       addition <- matrix(1, nrow=4, ncol=increment)
       iparents <- cbind(iparents, addition)
     }

     iparents_r <- dim(iparents)[2]
     dim_aij <- dim(listA)[1]
     rcA <- t(listA[, 1:2])
     aij <- listA[, 3]/2
     answer <- as.single(rep(0,iparents_r))
       rm("listA")
     cat(paste("starting to make D..."))

     For_out <- .Fortran(dij, iparents=as.integer(iparents), iparents_r=as.integer(iparents_r), dim_aij=as.integer(dim_aij), rcA=as.integer(rcA), aij=as.single(aij), N=as.integer(N), answer=as.single(answer))$answer
     rm("rcA", "aij", "iparents")

     if(!increment==0){
       Dijs <- For_out[-c((length(For_out)-increment+1):length(For_out))]
       } else{Dijs <- For_out}
     rm("For_out")
   } else{
        require(multicore)

        cat(paste("starting to make D..."))
        Dijs <- pvec(seq(1,dim(iparents)[2],1), FUN = wrap_dij, parents = iparents, relatedness = listA, N = N, mc.set.seed = FALSE, mc.silent = FALSE, mc.cores = ncores, mc.cleanup = TRUE)
        rm("listA") 
      }
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


