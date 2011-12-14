DFC <- function(ij.grandparents){
   imom.gp <- ij.grandparents[c(1,2)]
   idad.gp <- ij.grandparents[c(3,4)]
   jmom.gp <- ij.grandparents[c(5,6)]
   jdad.gp <- ij.grandparents[c(7,8)]

   if(all(identical(imom.gp, jmom.gp) & identical(idad.gp, jdad.gp)) | all(identical(imom.gp, jdad.gp) & identical(idad.gp, jmom.gp))) {
      return(1)
	} else {return(0)}
}

