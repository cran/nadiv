\name{makeS}
\alias{makeS}
\title{Creates the additive genetic relationship matrix for the sex chromosomes
}
\description{This is a preliminary function and may still contain some errors.  The function returns the inverse of the additive relationship matrix in sparse matrix format for the sex chromosomes (i.e., either X or Z).
}
\usage{
makeS(pedigree, heterogametic = "0", DosageComp = c(NULL, "ndc", "hoin"), returnS = FALSE)
}
\arguments{
  \item{pedigree }{
A pedigree where the columns are ordered ID, Dam, Sire, Sex}
  \item{heterogametic }{
Character indicating the label corresponding to the heterogametic sex used in the "Sex" column}
  \item{DosageComp }{
A character indicating whether a model of no dosage compensation ("ndc") or homozygous inactivation ("hoin") should be used.  The former ("ndc") follows the algorithm presented in Fernando and Grossman (1990) and has been validated.  However, the latter ("hoin") is still being validated and results should be used cautiously.  If "NULL" then the "ndc" model is assumed.
} 
  \item{returnS }{
Logical statement, indicating if the relationship matrix should be constructed in addition to the inverse}
}

\details{
Missing parents (e.g., base population) should be denoted by either 'NA', '0', or '*'.

The inverse of the sex-chromosome additive genetic relationship matrix (S-matrix) is constructed according to the algorithm presented in Fernando & Grossman (1990).  Additionally, the S-matrix itself can be constructed (although this takes much longer than computing S-inverse directly).

The algorithm in Fernando & Grossman (1990) does not account for any form of sex chromosome dosage compensation mechanism.  I have implemented one model, random X-inactivation (i.e., lyonization), in an attempt for more congruence with the biology of many organisms (e.g., eutherian mammals).  This second model of random inactivation in the homozygous sex is still preliminary and has yet to be properly peer reviewed. 
}


\value{
  \item{model }{the model of sex-chromosome dosage compensation assumed.}
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

