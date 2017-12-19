#' Combined p-value by minP and Fisher approached of aSPUw test incorporating multiple sources of functional annotations (aSPUmin)
#'
#' It gives p-values of the minP and Fisher approaches of multiple aSPUw tests.
#'
#' @param Y Response or phenotype data. It can be a quantitative trait. A vector with length n (number of observations).
#'
#' @param X Genotype or other data; each row for a subject, and each column
#'     for an SNP (or a predictor). The value of each SNP is the # of the copies
#'     for an allele. A matrix with dimension n by k (n : number of observation, k : number of SNPs (or predictors) ). No NAs are allowed. Please impute the NA vslues first.
#'
#' @param Z Covariates. A matrix with dimension n by p (n :number of observation, p : number of covariates).
#'
#' @param annot Functional annotations. A matrix with dimension k by r (k :number of SNPs, r : number of annotations).
#'
#' @param annot.wt The genome-wide informations. A vector with length r (number of annotations).
#'
#' @param pow SNP specific power(gamma values) used in FunSPU test.
#'
#' @param pow.annot Annotation specific power(gamma a values) used in FunSPU test.
#'
#' @param B number of permutations.
#'
#' @export
#' @return P-values for minP and Fisher approaches combining aSPUw tests.
#'
#' @author Yiding Ma and Peng Wei
#'
#' @references
#' Yiding Ma and Peng Wei (2017)
#' A Powerful and Adaptive Test for Incorporating Multiple Sources of Biological Knowledge into Association Analysis of Whole Genome Sequencing Data (Publication)
#'
#' @examples
#'
#' data(dat)
#' ## example analysis using aSPUwmin test on data.ex data.
#' test.aSPUwmin <- aSPUwmin(Y=data.ex$Y, X=data.ex$X, annot=data.ex$annot,
#'  annot.wt=data.ex$annot.wt, pow=c(1:6), B=1000)
#'
#' test.aSPUwmin$minP
#' # This is p-value of minP approach.
#'
#' test.aSPUwmin$Fisher
#' # This is p-value of Fisher approach.
#' @seealso \code{\link{aSPUmin}}


aSPUwmin <- function(Y, X, annot,annot.wt=NULL, cov = NULL,  pow=c(1:8, Inf), B=1000){

  n <- length(Y)
  if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)
  k <- ncol(X)

  m <- dim(annot)[2]

  if (is.null(cov)){
    ## NO nuisance parameters:
    XUs<-Xg <- X
    Xbar<-apply(Xg, 2, mean)
    subtract<-function(x, y) { x - y }
    Xgb=t(apply(Xg, 1, subtract, Xbar))

    r=Y-mean(Y)

    U<-as.vector( t(Xg) %*% r)

    #     if( model == "binomial" ) {
    #       CovS <- mean(Y)*(1-mean(Y))*(t(Xgb) %*% Xgb)
    #     } else {
    #       CovS <- var(Y)*(t(Xgb) %*% Xgb)
    #     }

    ##cov of the score stats:
    CovS <- var(Y)*(t(Xgb) %*% Xgb)

  } else {
    ## with nuisance parameters:
    tdat1<-data.frame(trait=Y, cov)
    fit1<-glm(trait~.,family=model,data=tdat1)
    pis<-fitted.values(fit1)
    Us<-matrix(0, nrow=n, ncol=k)
    for(i in 1:k){
      tdat2<-data.frame(X1=X[,i], cov)
      fit2<-glm(X1~.,data=tdat2)
      X1mus<-fitted.values(fit2)
      r <- Y - pis
      Us[, i]<-(Y - pis)*(X[,i] - X1mus)
    }
    U<-apply(Us, 2, sum)
    CovS<-matrix(0, nrow=k, ncol=k)
    for(i in 1:n)
      CovS<-CovS + Us[i,] %*% t(Us[i,])
  }

  Vs<-diag(CovS)
  diagSDs<-ifelse(Vs>1e-20, sqrt(Vs), 1e-10)

  Ts = matrix(NA ,nrow=length(pow), ncol=m)

  for (j in 1:length(pow)){
    for(k in 1:m){
      if (pow[j]<Inf) Ts[j,k] = sum((annot[,k]*U/diagSDs)^pow[j])/annot.wt[k] else Ts[j,k] = max(abs(annot[,k]*U/diagSDs))/annot.wt[k]
    }
  }

  pPerm0 = matrix(NA ,nrow=length(pow), ncol=m)
  T0s = matrix(NA ,nrow=B, ncol=m)
  P0s = matrix(NA ,nrow=B, ncol=m)
  minp0 = matrix(NA ,nrow=B, ncol=m)

  s <- sample(1:10^5,1)

  for (j in 1:length(pow)){
    set.seed(s) # to ensure the same samples are drawn for each pow
    for (b in 1:B){
      r0 <- sample(r, length(r))
      U0 <- as.vector(t(XUs) %*% r0)

      for ( k in 1:m ){
        if (pow[j] < Inf){ T0s[b,k] = sum((annot[,k]*U0/diagSDs)^pow[j])/annot.wt[k] }
        if (pow[j] == Inf) {T0s[b,k] =  max(abs(annot[,k]*U0/diagSDs))/annot.wt[k] }
      }
    }

    for ( k in 1:m ){
      pPerm0[j, k] = sum(abs(Ts[j, k])<=abs(T0s[,k])) / B
      P0s[,k] = ( (B-rank(abs(T0s[,k]))) + 1 )/(B)
      if (j==1) minp0[,k]=P0s[,k] else minp0[which(minp0[,k]>P0s[,k]), k]=P0s[which(minp0[,k]>P0s[,k]), k]
    }
  }

  minP.T <- apply(minp0, 1, min)
  fisher.fun <- function(x) -2*sum(log10(x))
  Fisher.T <- apply(minp0, 1, fisher.fun)

  p0 <- apply(pPerm0, 2, min) #pPerm0=

  minP.0 <- min(p0) #p0=pPerm0(in aSPU)
  Fisher.0 <- fisher.fun(p0)

  p.minP <- (sum(minP.0>=minP.T)+1)/(B+1)
  p.fisher <- (sum(Fisher.0<=Fisher.T)+1)/(B+1)

  list(minP = p.minP , Fisher = p.fisher)
}
