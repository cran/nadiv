\name{findDFC}
\alias{findDFC}
\title{
Finds the double first cousins in a pedigree}
\description{
Given a pedigree, all pairs of individuals that are double first cousins are returned.}

\usage{
findDFC(pedigree)
}

\arguments{
  \item{pedigree }{
A pedigree with columns organized: ID, Dam, Sire}
}

\value{
  \item{PedPositionList }{gives the list of row numbers for all the pairs of indidivuals that are related as double first cousins}
  \item{DFC }{gives the list of IDs, as characters, for all the pairs of individuals that are related as double first cousins}
  \item{FamilyCnt }{If two individuals, i and j, are double first cousins, then i's siblings will also be double first cousins with j's siblings.  Therefore, this is the total number of family pairs where offspring are related as double first cousins.}
}
\author{\email{matthewwolak@gmail.com}
}

