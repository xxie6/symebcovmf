#' Compute estimate for the residual variance
#' 
#' @param S An n-by-n Gram/covariance matrix
#' @param L An n-by-K matrix
#' @param lambda A K-length vector of weights
#' @param n The number of samples
#' @param K The number of factors
#' @param R2 The expectation of the squared Frobenius norm of the residuals
#' 
#' @return The residual variance
#' @export
#' 
estimate_resid_s2 <- function(S = NULL, L = NULL, lambda = NULL, n, K = NULL, R2 = NULL){
  if (is.null(R2) == TRUE){
    EvTEv <- apply(L, 2, function(x){return(sum(x^2))})
    obj_func_fit <- sum((S - tcrossprod(L %*% diag(sqrt(lambda), ncol = K)))^2) + 
      sum(lambda^2) - sum(EvTEv^2 * lambda^2)
    resid_s2 <- obj_func_fit/(n*(n+1))
  } else {
    resid_s2 <- R2/(n*(n+1))
  }
  return(resid_s2)
}