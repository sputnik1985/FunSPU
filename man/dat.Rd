\docType{data}
\name{dat}
\alias{dat}
\title{An example data set}
\usage{
data(dat)
}
\description{
  The dat data set is list of five objects. data.ex$Y is a vector of length
  1000. data.ex$X is a matrix of 1000 by 20. This
  simulated X matrix is assumed to be rare variants (MAF=0.01). All
  elements are 0, 1 or 2. data.ex$annot is a matrix of 20 by 7. data.ex$annot.wt is a vector of length
  7. data.ex$associndx is a vector of length 20. The
  values are 1 or 0, 0 indicate corresponding column of X matrix have no
  association with Y vector.
}
\examples{
data(dat)
data.ex$X[c(1:10, 501:510),]
data.ex$Y[c(1:10, 501:510)]
data.ex$annot
data.ex$annot.wt

# 11th to 20th column of X have no association with Y
data.ex$associndx

}
\keyword{datasets}
