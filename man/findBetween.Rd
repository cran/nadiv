\name{findBetween}
\alias{findBetween}
\title{Determines how close estimates are to a true value}

\description{
Function used in \code{proLik} when constructing profile likelihoods.  This determines between which variance estimates a confidence limit sits and if this is within a certain threshold.
}

\usage{findBetween(x, output, side = "L", threshb)}

\arguments{
  \item{x }{the critical value to search for}
  \item{output }{output from \code{constrainFun}}
  \item{side }{a character indicating whether the "L" lower or "U" upper confidence limits should be searched}
  \item{threshb }{the threshold for determining if an estimate is close enough to \code{x}}
}

\details{Used within \code{proLik} function}

\author{\email{matthewwolak@gmail.com}}

\seealso{See Also \code{\link{proLik}}}
