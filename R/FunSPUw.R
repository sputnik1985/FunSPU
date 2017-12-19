#' Inverse variance weighted adaptive SPU test incorporating multiple ources of functional annotations (FunSPUw)
#'
#' It gives p-values of the FunSPUw test.
#'
#' @param Y Response or phenotype data. It can be a quantitative trait. A vector with length n (number of observations).
#'
#' @param X Genotype or other data; each row for a subject, and each column
#'     for an SNP (or a predictor). The value of each SNP is the # of the copies
#'     for an allele. A matrix with dimension n by k (n : number of observation, k : number of SNPs (or predictors) ). No NAs are allowed. Please impute the NA vslues first.
#'
#' @param Z Covariates. A matrix with dimension n by p (n :number of observation, p : number of covariates).
#'
#' @param annot Functional annotations. AA matrix with dimension k by r (k :number of SNPs, r : number of annotations).
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
#' @return P-values for SPUw and FunSPUw test.
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
#' ## example analysis using FunSPUw test on data.ex data.
#' test.FunSPUw <- FunSPUw(Y=data.ex$Y, X=data.ex$X, annot=data.ex$annot,
#'  annot.wt=data.ex$annot.wt, pow=c(1:6), pow.annot=c(1,2,4,8,Inf), B=1000)
#'
#' test.FunSPUw$Ts
#' # This is a matrix of test Statistics for SPUw tests.
#' # pow = 1:6 and pow.annot = 1,2,4,8,Inf
#' # So, there are 6*5 = 30 SPU values.
#' # FunSPU(i,j) corresponds pow = i , pow.annot = j
#'
#' test.FunSPUw$pvs
#' # This is a vector of p-values for SPUw tests.
#'
#' test.FunSPUw$aTs
#' # This is minimum of them, FunSPUw statistic.
#'
#' test.FunSPUw$aSPU
#' # This is p-value of FunSPUw statistic.
#' @seealso \code{\link{FunSPU}}


FunSPUw <- function(Y, X, Z=NULL, annot, pow=c(1:6), annot.wt=NULL,pow.annot=c(1,2,4,8, Inf), B=1000){

  ##sample size
  n <- length(Y)
  if (is.null(X) && length(X)>0) X=as.matrix(X, ncol=1)

  ##number of SNPs
  k <- ncol(X)
  ##number of annotations
  m <- dim(annot)[2]

  if (is.null(Z)){
    ## NO nuisance parameters:
    Xg <- X
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
    tdat1<-data.frame(trait=Y, Z)
    fit1<-glm(trait~.,family=model,data=tdat1)
    pis<-fitted.values(fit1)
    Us<-matrix(0, nrow=n, ncol=k)
    for(i in 1:k){
      tdat2<-data.frame(X1=X[,i], Z)
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

  ############################
  ####    annotSPU part ######
  ############################

  sum.U.K <- NA
  Ts.annot<- rep(NA, m)

  Ts1=matrix(NA,length(pow)*dim(annot)[2])

  for(i in 1:length(pow)){
    for ( k in 1:m){
      if (pow[i]<Inf) {
        sum.U.k <- sum((annot[,k]*U/diagSDs)^pow[i])
        Ts.annot[k] <- sign(sum.U.k)*(abs(sum.U.k) )^(1/pow[i])/annot.wt[k]
        # Ts.annot[k] <- sum((annot[,k]*U/diagSDs)^pow[i])
      }
      else Ts.annot[k] = max(abs(annot[,k]*U/diagSDs))/annot.wt[k]
    }
    Ts1[((i-1)*m+1):(i*m)]=Ts.annot
  }

  ###################################
  ####    residual permutation ######
  ###################################

  sum.U0.K <- NA
  T0s1 = matrix(0, nrow=B, ncol=length(pow)*dim(annot)[2])
  #   s <- sample(1:10^5,1)

  for (b in 1:B){
    #     set.seed(s)
    r0 <- sample(r, length(r))
    U0 <- as.vector(t(Xg) %*% r0)

    for(i in 1:length(pow)){
      T0s.annot<- rep(NA, m)
      for ( k in 1:m){
        if (pow[i]<Inf) {
          sum.U0.k<- sum((annot[,k]*U0/diagSDs)^pow[i])
          T0s.annot[k] <- sign(sum.U0.k)*(abs(sum.U0.k) )^(1/pow[i])/annot.wt[k]
          # T0s.annot[k] <- sum((annot[,k]*U0/diagSDs)^pow[i])
        }
        else T0s.annot[k] = max(abs(annot[,k]*U0/diagSDs))/annot.wt[k]
      }
      T0s1[b,((i-1)*m+1):(i*m)]=T0s.annot
    }
  }

  Ts2<-rep(0, length(pow)*length(pow.annot))
  T0s2<-matrix(0, nrow=B, ncol=length(pow)*length(pow.annot))
  for(j2 in 1:length(pow.annot)){
    for(j in 1:length(pow)){
      if (pow.annot[j2]<Inf) {
        Ts2[(j2-1)*length(pow) +j] = sum(Ts1[((j-1)*m+1):(j*m)]^pow.annot[j2])
      }
      else Ts2[(j2-1)*length(pow) +j] = max(abs(Ts1[((j-1)*m+1):(j*m)]))

      for(b in 1:B){
        if (pow.annot[j2]<Inf) {
          T0s2[b, (j2-1)*length(pow) +j] = sum(T0s1[b, ((j-1)*m+1):(j*m)]^pow.annot[j2])
        }
        else T0s2[b, (j2-1)*length(pow) +j] = max(abs(T0s1[b, ((j-1)*m+1):(j*m)]))
      }
    }
  }

  Ts <- matrix(data=Ts2, nrow=length(pow), ncol=length(pow.annot), byrow=F )

  # permutation-based p-values:
  pPerm2 = rep(NA, length(pow)*length(pow.annot));
  pvs = NULL;

  P0s2 <- matrix(0, nrow=B, ncol=length(pow)*length(pow.annot))
  for(j in 1:(length(pow)*length(pow.annot)) ){
    pPerm2[j] = sum( abs(Ts2[j]) < abs(T0s2[,j]))/B
    P0s2[,j] = (B-rank(T0s2[,j])+1)/B
  }
  pvs <- matrix(data=pPerm2, nrow=length(pow), ncol=length(pow.annot), byrow=F )

  minP0s2 = apply(P0s2, 1, min)
  minP2 =  (sum( min(pPerm2) > minP0s2 )+1)/(B+1)
  cat("P0s caculated","\n")

  colnames(Ts) <- paste("SPU_out", pow.annot, sep="")
  rownames(Ts) <- paste("SPU_in", pow, sep="")

  colnames(pvs) <- paste("p_out", pow.annot, sep="")
  rownames(pvs) <- paste("p_in", pow, sep="")

  list(Ts = Ts, pvs = pvs, aTs=min(pPerm2), aSPU=minP2)
}
