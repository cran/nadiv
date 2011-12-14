\name{aiFun}
\alias{aiFun}
\title{Sampling (co)variances and correlations of random terms
}
\description{
This function returns the sampling covariances, variances, and correlations of the random effects fitted in an ASReml-R model}
\usage{
aiFun(model, Dimnames = NULL)
}
\arguments{
  \item{model }{
A model object returned by a call to the \code{asreml} function.}
  \item{Dimnames }{
A vector of characters if names are desired for the (co)variance/correlation matrix.}
}
\details{
The inverse of the average information provides the sampling (co)variance of each random term in the mixed model.  This function extracts the AI matrix from an ASReml-R model and organizes it so that the sampling covariances between random terms are located below the diagonal, the sampling variances of random terms are located along the diagonal, and the sampling correlations between random terms are located above the diagonal.  The order of the variances along the diagonal is the same as the order entered in the random section of the \code{asreml} function.  This is the same order as the rows of a call to the summary function, \code{summary(model)$varcomp}.}
\value{
A matrix of k x k dimensions is returned, if k is the number of random effects estimated in the model.  If \code{Dimnames} is specified, the row and column names are assigned according the vector of names in the argument.}
\references{
Gilmour, A.R., Gogel, B.J., Cullis, B.R., & Thompson, R. 2009. ASReml User Guide Release 3.0. VSN International Ltd., Hemel Hempstead, UK.}

\author{\email{matthewwolak@gmail.com}}

\note{
The vector of \code{Dimnames} should match the same order of variance components specified in the model (the row order in the \code{summary(model)$varcomp} object).}

\examples{
  \dontrun{
    library(asreml)
    ginvA <- asreml.Ainverse(warcolak)$ginv
    ginvD <- makeD(warcolak[,1:3])$listDinv
    warcolak$IDD <- warcolak$ID
    warcolak.mod <- asreml(phenotype ~ 1, random = ~ped(ID) + giv(IDD), ginverse = list(ID = ginvA, IDD = ginvD), data = warcolak) 
    summary(warcolak.mod)$varcomp
    aiFun(warcolak.mod, Dimnames = c("Va", "Vd", "Ve"))    
   }
}
