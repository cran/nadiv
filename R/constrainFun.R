constrainFun <- function(parameter.val, full, fm2, comp)
  {
  row <- which(fm2$Gamma == comp)
  fm2[row, 2:3] <- c(parameter.val, "F")
  full$G.param <- fm2
  con.mod <- update.asreml(object = full, trace = FALSE)
  if(con.mod$converge) return((-2*(full$loglik - con.mod$loglik))) else return(-998)
 }

