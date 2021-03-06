% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FunSPU.R
\name{FunSPU}
\alias{FunSPU}
\title{Adaptive SPU test incorporating multiple sources of functional annotations (FunSPU)}
\usage{
FunSPU(Y, X, Z = NULL, annot, annot.wt = NULL, pow = c(1:6),
  pow.annot = c(1, 2, 4, 8, Inf), B = 1000)
}
\arguments{
\item{Y}{Response or phenotype data. It can be a quantitative trait. A vector with length n (number of observations).}

\item{X}{Genotype or other data; each row for a subject, and each column
for an SNP (or a predictor). The value of each SNP is the # of the copies
for an allele. A matrix with dimension n by k (n : number of observation, k : number of SNPs (or predictors) ). No NAs are allowed. Please impute the NA vslues first.}

\item{Z}{Covariates. A matrix with dimension n by p (n :number of observation, p : number of covariates).}

\item{annot}{Functional annotations. A matrix with dimension k by r (k :number of SNPs, r : number of annotations).}

\item{annot.wt}{The genome-wide informations. A vector with length r (number of annotations).}

\item{pow}{SNP specific power(gamma values) used in FunSPU test.}

\item{pow.annot}{Annotation specific power(gamma a values) used in FunSPU test.}

\item{B}{number of permutations.}
}
\value{
P-values for SPU and FunSPU test.
}
\description{
It gives p-values of the FunSPU test.
}
\examples{

data(dat)
## example analysis using FunSPU test on data.ex data.
test.FunSPU <- FunSPU(Y=data.ex$Y, X=data.ex$X, annot=data.ex$annot,
 annot.wt=data.ex$annot.wt, pow=c(1:6), pow.annot=c(1,2,4,8,Inf), B=1000)

test.FunSPU$Ts
# This is a matrix of test Statistics for SPU tests.
# pow = 1:6 and pow.annot = 1,2,4,8,Inf
# So, there are 6*5 = 30 SPU values.
# FunSPU(i,j) corresponds pow = i , pow.annot = j

test.FunSPU$pvs
# This is a vector of p-values for SPU tests.

test.FunSPU$aTs
# This is minimum of them, FunSPU statistic.

test.FunSPU$aSPU
# This is p-value of FunSPU statistic.

}
\references{
Yiding Ma and Peng Wei (2017)
A Powerful and Adaptive Test for Incorporating Multiple Sources of Biological Knowledge into Association Analysis of Whole Genome Sequencing Data (Publication)
}
\seealso{
\code{\link{FunSPUw}}
}
\author{
Yiding Ma and Peng Wei
}
