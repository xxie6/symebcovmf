#' Compute the expectation of the squared Frobenius norm of the residuals
#' 
#' @param S An n-by-n Gram/covariance matrix
#' @param L An n-by-K matrix
#' @param lambda A K-length vector of weights
#' @param K The number of factors
#' 
#' @return The expectation of the squared Frobenius norm of the residuals
#' @export
#' 
compute_R2 <- function(S, L, lambda, K){
  EvTEv <- apply(L, 2, function(x){return(sum(x^2))})
  R2 <- sum((S - tcrossprod(L %*% diag(sqrt(lambda), ncol = K)))^2) + 
    sum(lambda^2) - sum(EvTEv^2 * lambda^2)
  return(R2)
}