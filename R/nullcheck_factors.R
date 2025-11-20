#' Check for factors that are zero vectors or have zero weight
#'
#' @param sym_ebcovmf_obj A symEBcovMF object
#' @param L2_tol The threshold for the L2 norm of the factor
#'
#' @return A symEBcovMF object
#' @export
#'
nullcheck_factors <- function(sym_ebcovmf_obj, L2_tol = 10^(-8)){
  null_lambda_idx <- which(sym_ebcovmf_obj$lambda == 0)
  factor_L2_norms <- apply(sym_ebcovmf_obj$L_pm, 2, function(v){sqrt(sum(v^2))})
  null_factor_idx <- which(factor_L2_norms < L2_tol)
  null_idx <- unique(c(null_lambda_idx, null_factor_idx))

  keep_idx <- setdiff(c(1:length(sym_ebcovmf_obj$lambda)), null_idx)

  if (length(keep_idx) < length(sym_ebcovmf_obj$lambda)){
    #remove factors
    sym_ebcovmf_obj$L_pm <- sym_ebcovmf_obj$L_pm[,keep_idx]
    sym_ebcovmf_obj$lambda <- sym_ebcovmf_obj$lambda[keep_idx]
    sym_ebcovmf_obj$KL <- sym_ebcovmf_obj$KL[keep_idx]
    sym_ebcovmf_obj$fitted_gs <- sym_ebcovmf_obj$fitted_gs[keep_idx]
  }

  #shouldn't need to recompute objective function or other things
  return(sym_ebcovmf_obj)
}
