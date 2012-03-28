makeA <- function(pedigree, method = NULL)
{

  numped <- numPed(pedigree)
  N <- dim(numped)[1]
  if(is.null(method)) method <- "fast"
  if(method == "fast"){
    numped[numped == -998] <- NA
    pedmm <- pedigree(sire = numped[,3], dam = numped[,2], label = numped[,1])
    tL <- relfactor(pedmm, pedmm@label)
    A <- as(crossprod(tL), "dgCMatrix")
    rm("tL")
  }


  if(method == "lowestmem"){
    ped <- t(numped)
    Tn <- (N*(N+1))/2 
    T <- rep(0, Tn)
    Dii <- rep(0, N)
  
    out <- .Fortran(a_ml, N=as.integer(N),ped=as.integer(ped),Tn=as.integer(Tn), T=as.single(T), Dii=as.single(Dii))

    nonzero.T <- which(out$T !=0)
    T.row <- as.integer(unlist(mapply(rep, 1:N, 1:N))[nonzero.T])
    T.col <- as.integer(unlist(mapply(seq, 1, length.out = 1:N))[nonzero.T])
    T.x <- out$T[nonzero.T]
    T <- sparseMatrix(i = as.integer(T.row), j = as.integer(T.col), x = T.x)
 
    D <- Diagonal(N, sqrt(out$Dii))
    rm("out")
    L <- suppressMessages(T%*%D)
    rm("T", "D")
    A <- as(tcrossprod(L), "dgCMatrix")
    rm("L")
  }  


  if(method == "flowmem"){
    numpedmm <- numped
    numpedmm[numped == -998] <- NA
    pedmm <- pedigree(sire = numpedmm[,3], dam = numpedmm[,2], label = numpedmm[,1])
    dnmiss <- which(numped[, 2] != -998)
    snmiss <- which(numped[, 3] != -998)
    Tinv.row <- c(numped[, 1][dnmiss], numped[, 1][snmiss], 1:N)
    Tinv.col <- c(numped[, 2][dnmiss], numped[, 3][snmiss], 1:N)
    Tinv.x <- c(rep(-0.5, length(dnmiss) + length(snmiss)), rep(1,N))
    el.order <- order(Tinv.col + Tinv.row/(N + 1), decreasing = FALSE)

    Tinv <- Matrix(0, N, N, sparse = TRUE)
    Tinv[1,2] <- 1
    Tinv@i <- as.integer(Tinv.row[el.order] - 1)
    Tinv@p <- as.integer(c(match(1:N, Tinv.col[el.order]), length(el.order) + 1) - 1)
    Tinv@x <- as.double(Tinv.x[el.order])
    dii <- Dmat(pedmm)
    D <- Diagonal(length(dii), sqrt(dii))
    rm("dii")
    Tt <- solve(Tinv)
    rm("Tinv")
    Tn <- dim(Tt)[1]
    Tlist <- data.frame(row = rep(1:Tn, Tn), column = rep(1:Tn, each = Tn), T = Tt@x)
    Tlist <- Tlist[as.vector(lower.tri(Tt, TRUE)),]
    Tlist <- Tlist[which(Tlist$T != 0),]
    order.index<-order(Tlist$column + Tlist$row/(Tn+1), decreasing=FALSE)
    T <-new("dtTMatrix", uplo = "L", i = as.integer(Tlist$row-1), j = as.integer(Tlist$column-1), x = Tlist$T, Dim = c(Tn, Tn))
    rm("Tlist")
    L <- suppressMessages(T%*%D)
    rm("T", "D")
    A <- as(tcrossprod(L), "dgCMatrix")
    rm("L")
  }

A
}

