#' Fit a rank K symEBcovMF model by greedily adding factors
#'
#' @param sym_ebcovmf_obj A symEBcovMF object
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
sym_ebcovmf_greedy <- function(sym_ebcovmf_obj, ebnm_fn, Kmax, maxiter, rank_one_tol, tol, sign_constraint = 'nonnegative', refit_lam = FALSE){
  # check if the symebcovmf object has a non-null L estimate
  if (is.null(sym_ebcovmf_obj$L_pm) == FALSE){
    R <- sym_ebcovmf_obj$S - tcrossprod(sym_ebcovmf_obj$L_pm %*% diag(sqrt(sym_ebcovmf_obj$lambda), ncol = ncol(sym_ebcovmf_obj$L_pm)))
  }
  
  # initialize greedy fitting procedure
  curr_rank <- 0
  obj_diff <- Inf
  while ((curr_rank < Kmax) & (obj_diff > tol)){
    # add factor
    print(paste('Adding factor', (curr_rank + 1)))
    sym_ebcovmf_obj <- sym_ebcovmf_r1_fit(R, 
                                          sym_ebcovmf_obj, 
                                          ebnm_fn, 
                                          maxiter, 
                                          rank_one_tol, 
                                          sign_constraint = sign_constraint)

    # check if new factor was added; if not then break
    if (length(sym_ebcovmf_obj$vec_elbo_K) == curr_rank){
      print(paste('Adding factor', (curr_rank + 1), 'does not improve the objective function'))
      break
    } else {
      # refit lambda if refit_lam = TRUE
      if (curr_rank > 0){
        if (refit_lam == TRUE){
          sym_ebcovmf_obj <- refit_lambda(sym_ebcovmf_obj)
        }
        # compute the objective function difference
        obj_diff <- sym_ebcovmf_obj$vec_elbo_K[curr_rank + 1] - sym_ebcovmf_obj$vec_elbo_K[curr_rank]
        if (obj_diff < tol){
          print(paste0('The objective function difference is: ', obj_diff, '. We will stop adding factors.' ))
        }
      }
    }
    curr_rank <- curr_rank + 1
  }

  return(sym_ebcovmf_obj)
}