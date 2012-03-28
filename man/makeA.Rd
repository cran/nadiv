\name{makeA}
\alias{makeA}
\title{Creates the additive genetic relationship matrix
}
\description{This returns the additive relationship matrix in sparse matrix format.
}
\usage{
makeA(pedigree, method = NULL)
}
\arguments{
  \item{pedigree }{
A pedigree where the columns are ordered ID, Dam, Sire}
  \item{method }{
Indication of method to use.  Either NULL, in which case the "fast" method is used by default, or one of "fast", "lowestmem", or "flowmem".  Method "fast" uses a good deal of RAM and is fastest for smaller pedigrees.  Method "lowestmem" is the slowest method, but uses the least amount of RAM.  Method "flowmem" is slower than "fast" for small pedigrees, but faster then "fast" for very large pedigrees (testing by the author puts the swtiching point at around 3000).  However, "flowmem" still uses more RAM than "lowestmem" for very large pedigress (despite being much faster).}
}
\details{
Missing parents (e.g., base population) should be denoted by either 'NA', '0', or '*'.

Uses \code{\link{relfactor}} function in package \code{pedigreemm}.  Also, used as a support function to \code{\link{makeD}}.}

\value{Returns A, or the numerator relationship matrix, in sparse matrix form.}

\author{\email{matthewwolak@gmail.com}}

\seealso{\code{\link{relfactor}}, \code{\link{makeD}}}

\examples{
 makeA(Mrode2)
}

