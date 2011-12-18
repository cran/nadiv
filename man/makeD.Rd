\name{makeD}
\alias{makeD}
\title{Creates the dominance genetic realationship matrix
}
\description{Given a pedigree, the matrix of coefficients of fraternity are returned - the D matrix.  Note, no inbreeding must be assumed.  Will return the inverse of the D matrix by default, otherwise this operation can be skipped if desired.
}
\usage{
makeD(pedigree, invertD = TRUE)
}
\arguments{
  \item{pedigree }{
A pedigree with columns organized: ID, Dam, Sire}

  \item{invertD }{
A logical indicating whether or not to invert the D matrix}
}
\details{
There exists no convenient method of obtaining the inverse of the dominance genetic relatedness matrix (or the D matrix itself) directly from a pedigree (such as for the inverse of A, i.e., Quaas (1995)).  Therefore, this function computes the coefficient of fraternity (Lynch and Walsh, 1998) for every individual in the pedigree with a non-zero additive genetic relatedness.  Note, the construction of the D matrix is more computationally demanding (in time and space) than is the construction of A.}
\value{
  \item{A }{the A matrix in sparse matrix form}
  \item{D }{the D matrix in sparse matrix form}
  \item{logDet }{the log determinant of the D matrix}
  \item{Dinv }{the inverse of the D matrix in sparse matrix form}
  \item{listDinv }{the three column form of the non-zero elements for the inverse of the D matrix}
}
\references{
Quaas, R.L. 1995. Fx algorithms. An unpublished note.

Lynch M., & Walsh, B. 1998. Genetics and Analysis of Quantitative Traits. Sinauer, Sunderland, Massachusetts.}
\author{\email{matthewwolak@gmail.com}
}

\seealso{
 \code{\link{makeDsim}} 
}


\examples{
  DinvMat <- makeD(Mrode9)$Dinv
}

