\name{makeDsim}
\alias{makeDsim}
\title{
Creates the dominance genetic relationship matrix through simulation}
\description{
Given a pedigree, the matrix of coefficients of fraternity are returned - the D matrix - as well as a simulated set of coefficients of fraternity (Ovaskainen et al. 2008).}
\usage{
makeDsim(pedigree, N, invertD = FALSE, calcSE = FALSE)
}

\arguments{
  \item{pedigree }{
A pedigree with columns organized: ID, Dam, Sire}
  \item{N }{
The number of times to simulate genotypes for the pedigree}
  \item{invertD }{
A logical indicating whether or not to invert the D matrix}
  \item{calcSE }{
A logical indicating whether or not the standard errors for each coefficient of fraternity should be calculated}
}

\details{
Missing parents (e.g., base population) should be denoted by either 'NA' or '0'.


Ovaskainen et al. (2008) indicated that the method of calculating the D matrix (see \code{\link{makeD}}) is only an approximation.  They proposed a simulation method that is implemented here.  This should be more appropriate when inbreeding occurs in the pedigree. 

The value, \code{listDsim} will list both the approximate values (returned from \code{\link{makeD}}) as well as the simulated values.  If \code{calcSE} is TRUE, these values will be listed in \code{listDsim}.
}

\value{
  \item{A }{the A matrix in sparse matrix form}
  \item{D }{the approximate D matrix in sparse matrix form}
  \item{logDetD }{the log determinant of the approximate D matrix}
  \item{Dinv }{the inverse of the approximate D matrix in sparse matrix form}
  \item{listDinv }{the three column form of the non-zero elements for the inverse of the approximate D matrix}
  \item{Dsim }{the simulated D matrix in sparse matrix form}
  \item{logDetDsim }{the log determinant of the simulated D matrix}
  \item{Dsiminv }{the inverse of the simulated D matrix in sparse matrix form}
  \item{listDsim }{the three column form of the non-zero and non-self elements for the simulated D matrix}
  \item{listDsiminv }{the three column form of the non-zero elements for the inverse of the simulated D matrix}

}
\references{
Ovaskainen, O., Cano, J.M., & Merila, J. 2008. A Bayesian framework for comparative quantitative genetics. Proceedings of the Royal Society B 275, 669-678.
}
\author{\email{matthewwolak@gmail.com}
}

\note{
This simulation can take a long time for large values of \code{N}.  If unsure, it is advisable to start with a lower N and gradually increase to gain a sense of the time required to execute a desired \code{N}.}

\seealso{
 \code{\link{makeD}} 
}


\examples{
  simDinvMat <- makeDsim(Mrode9, N = 1000, invertD = FALSE, calcSE = TRUE)$listDsim
}

