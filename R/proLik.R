proLik <- function(full.model, component, G = TRUE, negative = FALSE, nsample.units = 4, nse = 2, alpha = 0.05, threshold = 0.001, parallel = FALSE, ncores = getOption("cores"))
{

  df.in <- 1
  bound.check <- summary(full.model)$varcomp$constraint 
    if(any(bound.check == "Boundary")){
      warning("Boundary parameter: confidence interval estimated with a mix of Chi-squared distributions (df = 0.5)")
      df.in <- 0.5      	
      }

  s2 <- full.model$sigma2
  gamma.ind <- which(names(full.model$gammas) == component)
  gamma.est <- full.model$gammas[gamma.ind]
  std.err <- sqrt(aiFun(model = full.model)[gamma.ind, gamma.ind])
  thresh <- abs(gamma.est * threshold)
    low.limit <- (gamma.est - std.err * nse)
      if(gamma.est > 0 & low.limit < 0 & negative == FALSE) low.limit <- 0 
    bit1 <- seq(low.limit, gamma.est, length = nsample.units * 2)
    bit2 <- seq(gamma.est + (gamma.est*0.005), (gamma.est + std.err * (nse * 2)), length = (nsample.units * 3)+1)
  gamma.vec <- c(bit1, bit2)
  full.mod2 <- update.asreml(object = full.model, start.values = TRUE)$gammas.table

  if(parallel){
    require(multicore)
    profile <- list(lambdas = pvec(v = seq(1,length(gamma.vec),1), FUN = parConstrainFun, parameters = gamma.vec, full = full.model, fm2 = full.mod2, comp = component, G = G, mc.set.seed = FALSE, mc.silent = FALSE, mc.cores = ncores, mc.cleanup = TRUE), var.estimates = gamma.vec)
    } else{
    profile <- list(lambdas = vapply(gamma.vec, FUN = constrainFun, FUN.VALUE = vector("numeric", length = 1), full = full.model, fm2 = full.mod2, comp = component, G = G), var.estimates = gamma.vec)
      }

  chi.val <- -0.5 * qchisq(alpha, df = df.in, lower.tail = FALSE)
  leng.L <- length(bit1)
  leng.U <- length(bit2)

  pindex <- which(profile$lambdas == -998)
  pprob <- FALSE
    if(length(pindex) > 0){ 
      pprob <- TRUE
      Lpindex <- which(pindex <= leng.L)
      Upindex <- which(pindex > leng.U)
      leng.L <- leng.L - length(Lpindex)
      leng.U <- leng.U - length(Upindex)
      ptmp <- sapply(profile, "[", -pindex)
      profile <- list(lambdas = ptmp[,1], var.estimates = ptmp[,2])
    }

  

  if(profile$lambdas[leng.L + leng.U] > chi.val){
    if(pprob) warning("Some parameter values attempted did not converge: Upper confidence limit may not be estimable") else warning(paste(1-alpha, "% Upper confidence limit extends beyond the current profile likelihood: consider increasing value of 'nse'", sep = "")) 
    Upper.limit <- NULL

  } else{
      Uindex <- c((leng.L + 1):(leng.L + leng.U)) 
      Uprofile <- lapply(profile, "[", i = Uindex)
      UCI <- findBetween(chi.val, Uprofile, side = "U", threshb = thresh)
        if(UCI$more){
          UCI.tmp <- list(more = TRUE, hilo = NULL)
          Uprof.tmp <- list(lambdas = Uprofile$lambdas[UCI$hilo], var.estimates = Uprofile$var.estimates[UCI$hilo])
          while(UCI.tmp$more){
            tmp.gamma.vec <- seq(Uprof.tmp$var.estimates[1], Uprof.tmp$var.estimates[2], length = 8)
            if(parallel){
              tmp.lambdas <- pvec(v = seq(1,length(tmp.gamma.vec),1), FUN = parConstrainFun, parameters = tmp.gamma.vec, full = full.model, fm2 = full.mod2, comp = component, G = G, mc.set.seed = FALSE, mc.silent = FALSE, mc.cores = ncores, mc.cleanup = TRUE)
            } else{
              tmp.lambdas <- vapply(tmp.gamma.vec, FUN = constrainFun, FUN.VALUE = vector("numeric", length = 1), full = full.model, fm2 = full.mod2, comp = component, G = G)
             }

            UCI.tmp <- findBetween(x = chi.val, output = list(lambdas = tmp.lambdas, var.estimates = tmp.gamma.vec), side = "U", threshb = thresh)
            Uprof.leng <- length(Uprofile$lambdas)
            Uprofile$lambdas <- c(Uprofile$lambdas, tmp.lambdas)
	    Uprofile$var.estimates <- c(Uprofile$var.estimates, tmp.gamma.vec)
            if(UCI.tmp$more){
              Uprof.tmp <- list(lambdas = Uprofile$lambdas[Uprof.leng + UCI.tmp$hilo], var.estimates = Uprofile$var.estimates[Uprof.leng + UCI.tmp$hilo])
            } else break
          }
        Upper.limit <- UCI.tmp$estimate
        Uprofile$lambdas <- sort(Uprofile$lambdas, decreasing = TRUE)
        Uprofile$var.estimates <- sort(Uprofile$var.estimates, decreasing = FALSE)
        } else Upper.limit <- UCI$estimate
    }





  if(profile$lambdas[1] > chi.val){
    if(pprob) warning("Some parameter values attempted did not converge: Lower confidence limit may not be estimable") else warning(paste(1-alpha, "% Lower confidence limit extends beyond the current profile likelihood: consider increasing value of 'nse'", sep = "")) 
    Lower.limit <- NULL

  } else{
      Lindex <- c(1:(leng.L))
      Lprofile <- lapply(profile, "[", i = Lindex)
      LCI <- findBetween(chi.val, Lprofile, side = "L", threshb = thresh)
        if(LCI$more){
          LCI.tmp <- list(more = TRUE, hilo = NULL)
          Lprof.tmp <- list(lambdas = Lprofile$lambdas[LCI$hilo], var.estimates = Lprofile$var.estimates[LCI$hilo])
          while(LCI.tmp$more){
            tmp.gamma.vec <- seq(Lprof.tmp$var.estimates[2], Lprof.tmp$var.estimates[1], length = 8)
            if(parallel){
              tmp.lambdas <- pvec(v = seq(1,length(tmp.gamma.vec),1), FUN = parConstrainFun, parameters = tmp.gamma.vec, full = full.model, fm2 = full.mod2, comp = component, G = G, mc.set.seed = FALSE, mc.silent = FALSE, mc.cores = ncores, mc.cleanup = TRUE)
            } else{
                tmp.lambdas <- vapply(tmp.gamma.vec, FUN = constrainFun, FUN.VALUE = vector("numeric", length = 1), full = full.model, fm2 = full.mod2, comp = component, G = G)
                }
            LCI.tmp <- findBetween(x = chi.val, output = list(lambdas = tmp.lambdas, var.estimates = tmp.gamma.vec), side = "L", threshb = thresh)
            Lprof.leng <- length(Lprofile$lambdas)
            Lprofile$lambdas <- c(Lprofile$lambdas, tmp.lambdas)
	    Lprofile$var.estimates <- c(Lprofile$var.estimates, tmp.gamma.vec)
            if(LCI.tmp$more){
               Lprof.tmp <- list(lambdas = Lprofile$lambdas[Lprof.leng + LCI.tmp$hilo], var.estimates = Lprofile$var.estimates[Lprof.leng + LCI.tmp$hilo])
            } else break
          }
        Lower.limit <- LCI.tmp$estimate
        Lprofile$lambdas <- sort(Lprofile$lambdas, decreasing = FALSE)
        Lprofile$var.estimates <- sort(Lprofile$var.estimates, decreasing = FALSE)
        } else Lower.limit <- LCI$estimate
    }




  clambdas <- NULL
  cvar.estimates <- NULL
  if(is.null(Lower.limit)){
    if(any(pindex <= leng.L)){
      lo.seq <- c(1:leng.L)[-pindex[Lpindex]]
      clambdas <- profile$lambdas[lo.seq]
      cvar.estimates <- profile$var.estimates[lo.seq]
    }
  } else{
      clambdas <- Lprofile$lambdas
      cvar.estimates <- Lprofile$var.estimates
    }
  if(is.null(Upper.limit)){
    if(any(pindex > leng.L)){
      hi.seq <- c((leng.L + 1):(leng.L + leng.U))[-pindex[Upindex]]
      clambdas <- c(clambdas, profile$lambdas[hi.seq])
      cvar.estimates <- c(cvar.estimates, profile$var.estimates[hi.seq])
    }
  } else{
      clambdas <- c(clambdas, Uprofile$lambdas)
      cvar.estimates <- c(cvar.estimates, Uprofile$var.estimates)
    }


return(list(lambdas = clambdas, var.estimates = cvar.estimates * s2, UCL = Upper.limit * s2, LCL = Lower.limit * s2, component = component))
}

