\name{makeDomEpi}
\alias{makeDomEpi}
\title{Creates the additive by dominance and dominance by dominance epistatic genetic relationship matrices}
\description{
Given a pedigree, the matrix of additive by dominance (AD) genetic relatedness, dominance by dominance (DD) genetic relatedness, or both are returned.}
\usage{
makeDomEpi(pedigree, output = c("AD", "DD", "both"), Dinverse=FALSE)
}
\arguments{
  \item{pedigree }{
A pedigree where the columns are ordered ID, Dam, Sire}

  \item{output }{
Character(s) denoting which matrix and its inverse is to be constructed.}

  \item{Dinverse }{
A logical indicating whether or not to invert the D matrix}
}


\details{
Because of the computational demands of constructing the D matrix (see \code{\link{makeD}}), this function allows for the inverses that are derived from the D matrix (i.e., D-inverse, AD-inverse, and DD-inverse)to be constructed at the same time.  This way, the D matrix will only have to be constructed once for use in the three seperate genetic relatedness inverse matrices that depend upon it.  However, using the \code{output} and \code{Dinverse} options in different combinations will ensure that only the desired matrix inverses are constructed.

Both the AD and DD matrix are computed from the Hadamard product of the respective matrices (see also, \code{\link{makeAA}}).}

\value{All of the following will be returned.  However, the values of the \code{output} and \code{Dinverse} options passed to the function will determine which of the following are not NULL objects within the list:

  \item{D }{the D matrix in sparse matrix form}
  \item{logDetD }{the log determinant of the D matrix}
  \item{AD }{the AD matrix in sparse matrix form}
  \item{logDetAD }{the log determinant of the AD matrix}
  \item{DD }{the DD matrix in sparse matrix form}
  \item{logDetDD }{the log determinant of the DD matrix}
  \item{Dinv }{the inverse of the D matrix in sparse matrix form}
  \item{ADinv }{the inverse of the AD matrix in sparse matrix form}
  \item{DDinv }{the inverse of the DD matrix in sparse matrix form}
  \item{listDinv }{the three column form of the non-zero elements for the inverse of the D matrix}
   \item{listADinv }{the three column form of the non-zero elements for the inverse of the AD matrix}
   \item{listDDinv }{the three column form of the non-zero elements for the inverse of the DD matrix}
}

\author{\email{matthewwolak@gmail.com}
}

\seealso{\code{\link{makeA}}, \code{\link{makeD}}, \code{\link{makeAA}}
}


\examples{
  Boutput <- makeDomEpi(Mrode9, output = "b", Dinverse = FALSE)
  str(Boutput)
	
  DADoutput <- makeDomEpi(Mrode9, output = "AD", Dinverse = TRUE)
  str(DADoutput)
}

