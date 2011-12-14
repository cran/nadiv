FindDFC<-function(pedigree)
{
  numeric.pedigree <- numPed(pedigree)
  ped <- cbind(numeric.pedigree, genAssign(numeric.pedigree), rep(0, dim(numeric.pedigree)[1]))
  num.out <- ped[ped[,4] >= 2, ]
  ni <- dim(num.out)[1]
  maxid <- max(num.out[,1]) 

  i <- unlist(mapply(rep, (maxid - ni + 1):(maxid-1), each = seq((ni-1), 1)))
  j <- unlist(mapply(seq, (maxid - ni + 2):(maxid), maxid))

  ps <- cbind(i, j, numeric.pedigree[i, 2:3], numeric.pedigree[j, 2:3])
  ps.noFS <-ps[-which(ps[,3] == ps[,5] & ps[,4] == ps[,6]),]

  gps <- cbind(numeric.pedigree[ps.noFS[,3], 2:3], numeric.pedigree[ps.noFS[,4], 2:3], numeric.pedigree[ps.noFS[,5], 2:3], numeric.pedigree[ps.noFS[,6], 2:3]) 

  dfcs.vec<-apply(gps, MAR = 1, FUN = DFC)
  indexed <- cbind(ps.noFS, dfcs.vec)
  num.dfcs <- indexed[indexed[,7] == 1, 1:6]
  cnt.dfcs <- dim(unique(num.dfcs[, 3:6]))[1]
  dfcs.names<-cbind(as.character(pedigree[num.dfcs[,1], 1]), as.character(pedigree[num.dfcs[,2], 1]))

return(list(PedPositionList = num.dfcs[, 1:2], DFC = dfcs.names, FamilyCnt = cnt.dfcs))
}
