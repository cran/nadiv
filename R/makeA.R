makeA<-function (pedigree) 
{
  numeric.pedigree <- numPed(pedigree)
  N <- dim(numeric.pedigree)[1]
  ped <- t(numeric.pedigree)
  T <- matrix(0, N, N)
  F <- rep(0, N)
  Dii <- rep(0, N)

  out <- .Fortran(makea_ml, ped_n=as.integer(N),ped=as.integer(ped),T=as.single(T),F=as.single(F),Dii=as.single(Dii),Aout=as.single(T))


  T.i <- as.integer(rep(1:N, N)[which(out$T > 0)]-1)
  T.x <- out$T[which(out$T > 0)]
  T.j <- as.integer(rep(1:N, each=N)[which(out$T>0)]-1)
  T <- new("dtTMatrix", Dim=as.integer(c(N,N)))
  T@uplo <- "L"
  T@i <- T.i
  T@x <- T.x
  T@j <- T.j

  D <- Diagonal(N, sqrt(out$Dii))
  L <- T%*%D
  Asparse <- as(L%*%t(L), "dgCMatrix")
  Adense <- as(Asparse, "matrix")
  
return(list(Asparse = Asparse, Adense = Adense))
}

