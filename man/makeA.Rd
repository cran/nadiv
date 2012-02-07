\name{makeA}
\alias{makeA}
\title{Creates the additive genetic relationship matrix
}
\description{This returns the additive relationship matrix in sparse matrix format.
}
\usage{
makeA(pedigree)
}
\arguments{
  \item{pedigree }{
A pedigree where the columns are ordered ID, Dam, Sire}
}
\details{
Missing parents (e.g., base population) should be denoted by either 'NA' or '0'.

Uses \code{\link{relfactor}} function in package \code{pedigreemm}.  Also, used as a support function to \code{\link{makeD}}.}

\value{Returns A, or the numerator relationship matrix, in sparse matrix form.}

\author{\email{matthewwolak@gmail.com}}

\seealso{\code{\link{relfactor}}, \code{\link{makeD}}}

\examples{
 makeA(Mrode2)
}

