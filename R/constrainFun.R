constrainFun <- function(parameter.val, full, fm2, comp, G, mit = 600)
  {
  row <- which(fm2$Gamma == comp)
  fm2[row, 2:3] <- c(parameter.val, "F")
  if(G) full$G.param <- fm2 else full$R.param <- fm2
  con.mod <- update.asreml(object = full, maxiter = mit, trace = FALSE)
  if(con.mod$converge) return((-2*(full$loglik - con.mod$loglik))) else return(NA)
 }

