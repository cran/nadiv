\name{constrainFun}
\alias{constrainFun}
\title{
Function used in the \code{proLik} function to produce a profile likelihood for a variance component}

\description{
Given a model object from \code{asreml} and a range of estimates of the parameter, the function will supply the likelihood ratio test statistic for the comparison of the full model to one where the parameter of interest is constrained.}

\usage{
constrainFun(parameter.val, full, fm2, comp)
}

\arguments{
  \item{parameter.val }{a value for which the log Likelihood of a model is to be calculated}
  \item{full }{the full model \code{asreml} object}
  \item{fm2 }{starting values for the full model}
  \item{comp}{which variance component to constrain}
}

\details{
Used internally in the \code{proLik} function}

\author{\email{matthewwolak@gmail.com}}

\seealso{See Also \code{\link{proLik}}}
