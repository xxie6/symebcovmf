#' Compute the objective function
#' 
#' @param S An n-by-n Gram/covariance matrix
#' @param L An n-by-K matrix
#' @param lambda A K-length vector of weights
#' @param resid_s2 The residual variance
#' @param n The number of samples
#' @param K The number of factors
#' @param KL A K-length vector containing the negative KL divergence
#' @param R2 The expectation of the squared Frobenius norm of the residuals
#' 
#' @return The objective function
#' @export
#' 
compute_elbo <- function(S = NULL, L = NULL, lambda = NULL, resid_s2, n, K = NULL, KL, R2 = NULL){
  if (is.null(R2) == TRUE){
    EvTEv <- apply(L, 2, function(x){return(sum(x^2))})
    elbo <- -((n*(n-1)/4)*log(2*pi*resid_s2)) - 
      (n/2)*log(2*pi*2*resid_s2) -
      (1/(4*resid_s2))*(sum((S - tcrossprod(L %*% diag(sqrt(lambda), ncol = K)))^2) + 
                          sum(lambda^2) - sum(EvTEv^2 * lambda^2)) + sum(KL)
  } else {
    elbo <- -((n*(n-1)/4)*log(2*pi*resid_s2)) - (n/2)*log(2*pi*2*resid_s2) -
      (1/(4*resid_s2))*R2 + sum(KL)
  }
  return(elbo)
}