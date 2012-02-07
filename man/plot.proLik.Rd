\name{plot.proLik}
\alias{plot.proLik}
\alias{plot}

\title{Graphically depicts the profile Likelihood}

\description{This function plots the output from a \code{proLik} object.}

\usage{
  \method{plot}{proLik}(x, CL = TRUE, alpha = 0.05, type = "l", \dots)
}

\arguments{
  \item{x }{the output from \code{proLik}}
  \item{CL }{a logical indicating whether a line representing the Confidence Limit is to be drawn}
  \item{alpha }{if \code{CL} is TRUE, the confidence level at which \code{CL} is to be drawn}
  \item{type }{the type of plot to be generated, see arguments to \code{plot}}
  \item{\dots }{other arguments to \code{plot}}
}

\author{\email{matthewwolak@gmail.com}}


\seealso{\code{\link{proLik}}}

\examples{##See examples in \code{\link{proLik}}}
