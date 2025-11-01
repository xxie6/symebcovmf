#' Fit rank K symEBcovMF model
#'
#' @param S A n-by-n Gram/covariance matrix
#' @param ebnm_fn An EBNM solver
#' @param Kmax The maximum number of factors to add
#' @param maxiter The maximum number of iterations used when optimizing the rank one fit
#' @param rank_one_tol The convergence tolerance for the rank one fit
#' @param tol The convergence tolerance for the rank K fit
#' @param sign_constraint The sign constraint (if any) for the loadings estimate
#' @param refit_lam True or False for if you want to refit the lambda values after each factor is added
#'
#' @return A symEBcovMF object
#' @export
#'
sym_ebcovmf_fit <- function(S, ebnm_fn, Kmax, maxiter, rank_one_tol, tol, sign_constraint = 'nonnegative', refit_lam = FALSE){
  #initialize object
  sym_ebcovmf_obj <- sym_ebcovmf_init(S)

  curr_rank <- 0
  obj_diff <- Inf
  while ((curr_rank < Kmax) & (obj_diff > tol)){
    # add factor
    sym_ebcovmf_obj <- sym_ebcovmf_r1_fit(S, sym_ebcovmf_obj, ebnm_fn, maxiter, rank_one_tol, sign_constraint = sign_constraint)

    # check if new factor was added
    if (length(sym_ebcovmf_obj$vec_elbo_K) == curr_rank){
      print(paste('Adding factor', (curr_rank + 1), 'does not improve the objective function'))
      break
    } else {
      if (curr_rank > 0){
        if (refit_lam == TRUE){
          sym_ebcovmf_obj <- refit_lambda(S, sym_ebcovmf_obj)
        }
        obj_diff <- sym_ebcovmf_obj$vec_elbo_K[curr_rank + 1] - sym_ebcovmf_obj$vec_elbo_K[curr_rank]
      }
    }
    curr_rank <- curr_rank + 1
  }

  return(sym_ebcovmf_obj)
}
