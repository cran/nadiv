\name{warcolak}
\alias{warcolak}
\docType{data}
\title{
Pedigree and phenotypic values for a mythical half-sib population of Warcolaks}
\description{
An example pedigree with phenotypic data constructed from simulated additive, dominance, and environmental effects.}
\usage{data(warcolak)}
\format{
  A data frame with 480 observations on the following 4 variables.
  \describe{
    \item{\code{ID }}{a factor specifying 480 unique individual IDs}
    \item{\code{Dam }}{a factor specifying the unique ID for each of 3 dams per sire}
    \item{\code{Sire }}{a factor specifying the unique ID for all 30 sires}
    \item{\code{phenotype} }{a numeric vector of phenotypic values}
  }
}
\details{
The dataset was simulated to have a population mean of 50 (with no other fixed effects) and additive, dominance, and environmental variances of 40, 30, & 30, respectively.  The random effects were drawn from a multivariate random normal distribution [e.g., additive effects: N ~ (0, A * Va)] with means of zero and variances equal to the product of the desired additive genetic variance and the additive genetic relatedness matrix.  Because of this, the actual variance in random effects will vary slightly from the amount specified in the simulation (i.e., 40, 30, & 30) 
}

\keyword{datasets}
