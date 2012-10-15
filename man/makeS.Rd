\name{makeS}
\alias{makeS}
\title{Creates the additive genetic relationship matrix for the sex chromosomes
}
\description{First, this is a preliminary function and may still contain some errors.  The function returns the inverse of the additive relationship matrix in sparse matrix format for the sex chromosomes (i.e., either X or Z).
}
\usage{
makeS(pedigree, heterogametic = "0", returnS = FALSE)
}
\arguments{
  \item{pedigree }{
A pedigree where the columns are ordered ID, Dam, Sire, Sex}
  \item{heterogametic }{
Character indicating the label corresponding to the heterogametic sex used in the "Sex" column}
  \item{returnS }{
Logical statement, indicating if the relationship matrix should be constructed in addition to the inverse}
}

\details{
This is the first version of this function and so any problems should be reported to the author.

Missing parents (e.g., base population) should be denoted by either 'NA', '0', or '*'.

The inverse of the sex-chromosome additive genetic relationship matrix (S-matrix) is constructed according to the algorithm presented in Fernando & Grossman (1990).  Additionally, the S-matrix itself can be constructed (although this takes much longer than computing S-inverse directly). 
}


\value{
  \item{S }{the sex-chromosome relationship matrix in sparse matrix form or NULL if \code{returnS} = FALSE}
  \item{Sinv }{the inverse of the S matrix in sparse matrix form}
  \item{listSinv }{the three column form of the non-zero elements for the inverse of the S matrix}
}

\references{
Fernando, R.L. & Grossman, M. 1990. Genetic evaluation with autosomal and X-chromosomal inheritance. Theoretical and Applied Genetics, 80:75-80.
}
 

\author{\email{matthewwolak@gmail.com}}

\examples{
 makeS(FG90, heterogametic = "0", returnS = TRUE)
}

