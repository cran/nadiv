findBetween <- function(x, output, side = "L", threshb)
{
  if(side == "L"){
    lo <- max(which(output$lambdas < x)) 
    hi <- min(which(output$lambdas > x))
    } else{ 
        lo <- min(which(output$lambdas < x)) 
        hi <- max(which(output$lambdas > x))
      }
  if(!is.integer(hi)){
    stop(paste("Boundary of", if(side == "L") "Lower" else "Upper", "Confidence Limit cannot be determined: increase 'nsample.units' & maybe also consider decreasing 'nse'", sep = " "))
  }
  if(!is.integer(lo)){
    stop(paste("Boundary of", if(side == "L") "Lower" else "Upper", "Confidence Limit cannot be determined: increase 'nsample.units' & increase 'nse'", sep = " "))
  }
  slope <- (output$lambdas[hi] - output$lambdas[lo]) /
(output$var.estimates[hi] - output$var.estimates[lo])
  interp <- output$var.estimates[lo] + ((x - output$lambdas[lo]) / slope)
  diff.hi <- abs(output$var.estimates[hi] - interp)
  diff.lo <- abs(output$var.estimates[lo] - interp) 
  if(diff.hi <= threshb & diff.lo <= threshb){
    return(list(more = FALSE, estimate = interp))
    } else{
        hilo <- c(output$var.estimates[hi], output$var.estimates[lo])
        use <- c("hi", "lo")
        return(list(more = TRUE, hilo = c(hi, lo)))
      }
}

