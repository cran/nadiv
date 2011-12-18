\name{DFC}
\alias{DFC}
\title{Indicates Double First Cousins}
\description{
Given all eight grandparents of two individuals, the function will indicate whether they are double first cousins.  Mainly to be used inside of \code{\link{FindDFC}}}

\usage{
DFC(ij.grandparents)
}

\arguments{
  \item{ij.grandparents }{
a vector indicating the dam and sire of both i's and j's dam and sire}
}

\value{
returns either a 0 if individual's i and j are not double first cousins, or a 1 if they are}

\author{\email{matthewwolak@gmail.com}
}

\seealso{
 \code{\link{FindDFC}}
}

