\name{In.house}
\alias{In.house}
\docType{data}
\title{ IN-HOUSE NUCLEIC ACID AMPLIFICATION TESTS (INH) FOR TB PLEURITIS }


\description{
This dataset gives the observed cross-tabulation of the test under evaluation (nucleic acid amplification test)  and the reference test.  In studies 1 and 2, the reference test was a culture test, in studies 3 and 4 it was a composite reference test of culture and clinical data (including signs, symptoms and clinical response to empiric TB therapy), in study 5 it was a composite reference test of culture and microscopy and in studies 6 to 12 it was a composite reference test of culture and biopsy.
}

\usage{data(In.house)}


\format{
  A matrix with 12 observations on the following 4 variables.
  \describe{
    \item{\code{TP}}{Observed individuals who tested positive on both tests }
    \item{\code{FP}}{Observed individuals who tested positive on the test under evaluation and negative on the reference test}
    \item{\code{FN}}{Observed individuals who tested negative on the test under evaluation and positive on the reference test}
    \item{\code{TN}}{Observed individuals who tested negative on both tests }
  }
}


\references{
Pai, M. et al. (2004) \emph{Nucleic acid amplification in the diagnosis of tuberculous pleuritis: a systematic review and meta-analysis}. 
BMC Infect Dis 2004, 4:6.

}

\examples{
data(In.house)
In.house

}

\keyword{datasets}
