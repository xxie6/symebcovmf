#' Initialize the factors at specified values
#'
#' @param S A n-by-n Gram/covariance matrix
#' @param sym_ebcovmf_obj A symEBcovMF object
#' @param init_L An initial estimate for the loadings matrix. The columns of this matrix should have L2 norm equal to 1.
#' @param init_lambda A vector containing the initial values for the non-negative weights
#'
#' @return A symEBcovMF object
#' @export
#'
sym_ebcovmf_factors_init <- function(S, sym_ebcovmf_obj, init_L, init_lambda){
  init_L_norms <- apply(init_L, 2, function(x){sqrt(sum(x^2))})
  if(any(abs(init_L_norms - 1) > 10^(-8))){
    stop('The columns of the loadings initialization should have L2 norm equal to 1')
  }

  sym_ebcovmf_obj$L_pm <- init_L
  sym_ebcovmf_obj$lambda <- init_lambda

  sym_ebcovmf_obj$resid_s2 <- estimate_resid_s2(S = S,
                                                L = sym_ebcovmf_obj$L_pm,
                                                lambda = sym_ebcovmf_obj$lambda,
                                                n = nrow(init_L),
                                                K = length(sym_ebcovmf_obj$lambda))

  sym_ebcovmf_obj$elbo <- compute_elbo(S = S,
                                       L = sym_ebcovmf_obj$L_pm,
                                       lambda = sym_ebcovmf_obj$lambda,
                                       resid_s2 = sym_ebcovmf_obj$resid_s2,
                                       n = nrow(init_L),
                                       K = length(sym_ebcovmf_obj$lambda),
                                       KL = rep(0, length(sym_ebcovmf_obj$lambda)))
  return(sym_ebcovmf_obj)
}
