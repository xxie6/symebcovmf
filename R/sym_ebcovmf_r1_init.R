#' Initialize the greedily added factor
#'
#' @param R An n-by-n symmetric matrix
#' @param nonnegative A TRUE/FALSE value for whether the initialization is constrained to be nonnegative
#'
#' @return An initial estimate for the factor
#' @export
#'
sym_ebcovmf_r1_init <- function(R, nonnegative = TRUE){
  svd1 <- RSpectra::eigs_sym(R, k = 1)
  v <- svd1$vectors # scaled such that v'v = 1
  lambda_k <- svd1$values
  if(nonnegative == TRUE){
    svd_v <- v
    v <- pmax(svd_v, 0)
    if (sqrt(sum(v^2)) > 0){
      v <- v/sqrt(sum(v^2))
    }

    minus_v <- pmax(-svd_v, 0)
    if(sqrt(sum(minus_v^2)) > 0){
      minus_v <- minus_v/sqrt(sum(minus_v^2))
    }
  } else {
    minus_v <- -v
  }
  lambda_options <- c(drop(t(v) %*% R %*% v), drop(t(minus_v) %*% R %*% minus_v))
  if(which.max(lambda_options) == 1){
    v <- v
    lambda_k <- lambda_options[1]
  } else {
    v <- minus_v
    lambda_k <- lambda_options[2]
  }
  return(list(v = v, lambda_k = lambda_k))
}
