\name{proLik}
\alias{proLik}

\title{
Estimates the profile likelihood of a random effect
}

\description{
When a mixed model is run in ASReml-R, this function can estimate the profile likelihood of variance components from the output.  
}

\usage{
proLik(full.model, component, negative = FALSE, nsample.units = 4, nse = 2, alpha = 0.05, threshold = 0.05)
}

\arguments{
  \item{full.model }{An \code{asreml} model object}
  \item{component }{A character vector indicating for which variance component the profile likelihood will be constructed.  Must be an object in \code{full.model$gammas}.}
  \item{negative }{Logical indicating whether or not the \code{component} can take on a negative value (i.e., a covariance)}
  \item{nsample.units }{Number of sample units to be used in constructing the profile likelihood}
  \item{nse }{Number of standard errors on either side of the estimate, over which the profile likelihood should be constructed}
  \item{alpha }{The critical value for determining the Confidence Interval}
  \item{threshold }{Acceptable distance, between actual sample values and interpolated values, for determining the upper and lower limits of the Confidence Interval.  Actual samples will be no more than this distance from the true value of the estimate}
}


\details{
For the \code{negative} argument, this should be used if the profile likelihood of a covariance component is to be constructed.

The function uses a grid search to obtain the approximate confidence limits.  Therefore, \code{nsample.units} should be carefully thought about beforehand when running the function.  Increasing this value will ensure a smoother surface, but at cost to time.  The function weights this quantity to ensure more samples for regions of the profile likelihood greater than the estimate.  The lower bound of the sampling interval to the estimate is made up of 2 times \code{nsample.units}, whereas the sampling interval from the estimate to the upper bound is 3 times \code{nsample.units}.  

Similarly for \code{nse}, the function will include twice as many standard errors from the estimate to the upper bound than it does from the estimate to the lower bound.  If \code{negative} is FALSE, and the lower bound of the sampling interval extends beyond zero, this will instead be set to 0.
}


\value{
  \item{lambdas }{negative log Likelihood ratio test statistic.  Estimated from the log Likelihood of the \code{full.model} and the log Likelihood of the model with the \code{component} constrained to a value in the sampling interval}
  \item{var.estimates }{value along the sampling interval for which the \code{component} was constrained}
  \item{UCL }{approximate Upper Confidence Limit}
  \item{LCL }{approximate Lower Confidence Limit}
  \item{component }{the component for which the profile likelihood surface has been constructed}
}

\author{\email{matthewwolak@gmail.com}}

\section{Warning }{
May be unfeasible to estimate profile likelihoods for complex models with many variance components
}

\seealso{
\code{\link{aiFun}}
}

\examples{
  \dontrun{
    library(asreml)
    ginvA <- asreml.Ainverse(warcolak[, c(1,3,2)])$ginv
    ginvD <- makeD(warcolak[,1:3])$listDinv
    warcolak$IDD <- warcolak$ID
    warcolak.mod <- asreml(trait1 ~ 1, random = ~ped(ID) + giv(IDD), ginverse = list(ID = ginvA, IDD = ginvD), data = warcolak) 
    summary(warcolak.mod)$varcomp
    profileA <- proLik(full.model = warcolak.mod, component = "ped(ID)!ped", negative = FALSE, nsample.units = 3, nse = 3)
    profileA
    profileD <- proLik(full.model = warcolak.mod, component = "giv(IDD).giv", negative = FALSE, nsample.units = 3, nse = 3)

    x11(w = 6, h = 8)
    par(mfrow = c(2,1))
      plot.proLik(profileA) 
      plot.proLik(profileD)
   }

}


