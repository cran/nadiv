makeA <- function(pedigree)
{
  numped <- numPed(pedigree)
  N <- dim(numped)[1]
  dnmiss <- which(numped[, 2] != -998)
  snmiss <- which(numped[, 3] != -998)
  Tinv.row <- c(numped[, 1][dnmiss], numped[, 1][snmiss], 1:N)
  Tinv.col <- c(numped[, 2][dnmiss], numped[, 3][snmiss], 1:N)
  Tinv.x <- c(rep(-0.5, length(dnmiss) + length(snmiss)), rep(1, N))
  el.order <- order(Tinv.col + Tinv.row/(N + 1), decreasing = FALSE)
  bnmiss <- which(numped[, 2] != -998 & numped[, 3] != -998)
  selfed <- bnmiss[which((numped[bnmiss, 2] - numped[bnmiss, 3]) == 0)]
  if(length(selfed) > 0){
     repeat1 <- which(sapply(seq(2,length(Tinv.row)), FUN = function(x){identical(c(Tinv.row[el.order][x], Tinv.col[el.order][x]), c(Tinv.row[el.order][x-1], Tinv.col[el.order][x-1]))}))
     Tinv.row <- Tinv.row[el.order][-c(repeat1 + 1)]
     Tinv.col <- Tinv.col[el.order][-c(repeat1 + 1)]
     Tinv.x <- Tinv.x[el.order]
     Tinv.x[repeat1] <- -1
     Tinv.x <- Tinv.x[-c(repeat1 + 1)]
     Tinv <- Matrix(0, N, N, sparse = TRUE)
     Tinv[1, 2] <- 1
     Tinv@i <- as.integer(Tinv.row - 1)
     Tinv@p <- as.integer(c(match(1:N, Tinv.col), length(Tinv.col) + 1) - 1)
     Tinv@x <- as.double(Tinv.x)
  } else{
  Tinv <- Matrix(0, N, N, sparse = TRUE)
  Tinv[1, 2] <- 1
  Tinv@i <- as.integer(Tinv.row[el.order] - 1)
  Tinv@p <- as.integer(c(match(1:N, Tinv.col[el.order]), length(el.order) + 1) - 1)
  Tinv@x <- as.double(Tinv.x[el.order])
    }
  T <- as(solve(Tinv), "dgCMatrix")
  tT <- t(T)
  numped[numped == -998] <- N + 1
    Cout <- .C("ddiag",
	    as.integer(numped[, 2] - 1), 
	    as.integer(numped[, 3] - 1),
	    as.double(c(rep(0, N), -1)), 
	    as.double(rep(1, N)), 
            as.integer(tT@i), 
	    as.integer(c(tT@p, length(tT@x))), 
            as.double(tT@x), 
            as.integer(N), 
	    as.integer(length(tT@x)))
 as(T %*% Diagonal(N, Cout[[4]]) %*% tT, "symmetricMatrix")
}

