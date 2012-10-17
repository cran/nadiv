\name{numPed}
\alias{numPed}
\title{Creates a numeric form of a pedigree
}
\description{
This function takes a pedigree and converts it into a standard integer form used by many functions in the 'nadiv' package.}
\usage{
numPed(pedigree, check = TRUE)
}
\arguments{
  \item{pedigree }{
A three column pedigree object, where the columns correspond to: ID, Dam, & Sire}
  \item{check }{
A logical argument indicating if checks on the validity of the pedigree structure should be made}

}
\details{
Missing parents (e.g., base population) should be denoted by either 'NA', '0', or '*'.

Individuals must appear in the ID column in rows preceeding where they appear in either the Dam or Sire column.  If \code{check = FALSE} any pedigree will be transformed into a pedigree consisting of integers and missing values denoted by '-998'.  CAUTION: this should not be set to TRUE without careful thought.

Based on code from the \code{MCMCglmm} package}

\seealso{\code{\link[MCMCglmm]{MCMCglmm}}}


\value{
The pedigree, where individuals are now numbered from 1 to n and unknown parents are assigned a value of '-998'.
}


