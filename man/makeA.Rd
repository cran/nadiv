\name{makeA}
\alias{makeA}
\title{Creates the additive genetic relationship matrix
}
\description{This returns the additive relationship matrix in both sparse and dense matrix formats.
}
\usage{
makeA(pedigree)
}
\arguments{
  \item{pedigree }{
A pedigree where the columns are ordered ID, Dam, Sire}
}
\details{
Uses the algorithm of Meuwissen and Luo (1992).  Also, used as a support function to \code{\link{makeD}}.}

\value{Returns A, or the numerator relationship matrix, in sparse and dense matrix form.

  \item{Asparse }{A in sparse matrix form}
  \item{Adense }{A in dense matrix form}
}

\references{Meuwissen, T.H.E. & Luo, Z. 1992. Computing inbreeding coefficients in large populations. Genetics, Selection, Evolution 24, 305-313.}

\author{\email{matthewwolak@gmail.com}}

\seealso{\code{\link{makeD}}}

\examples{
 makeA(Mrode2)
}

