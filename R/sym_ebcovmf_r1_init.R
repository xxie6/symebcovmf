#' Initialize the greedily added factor
#' 
#' @param R An n-by-n symmetric matrix
#' 
#' @return An initial estimate for the factor
#' @export
#' 
sym_ebcovmf_r1_init <- function(R){
  svd1 <- RSpectra::eigs_sym(R, k = 1)
  v <- svd1$vectors # scaled such that v'v = 1
  if (abs(min(v)) > abs(max(v))){
    v <- -1*v
  }
  lambda_k <- svd1$values
  return(list(v = v, lambda_k = lambda_k))
}