#' Initialize rank K symEBcovMF Fit
#' 
#' @param S An n-by-n Gram/covariance matrix
#' 
#' @return A symEBcovMF object
#' @export
#' 
sym_ebcovmf_init <- function(S){
  sym_ebcovmf_obj <- list(n = ncol(S), 
                       L_pm = NULL, 
                       resid_s2 = 0, 
                       KL = c(), 
                       lambda = c(), 
                       elbo = 0, 
                       vec_elbo_K = c(), 
                       vec_elbo_full = c(), 
                       fitted_gs = list())
  return(sym_ebcovmf_obj)
}