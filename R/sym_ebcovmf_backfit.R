#' Fit a rank K symEBcovMF model using backfit
#'
#' @param S A n-by-n Gram/covariance matrix
#' @param sym_ebcovmf_obj A symEBcovMF object
#' @param ebnm_fn An EBNM solver
#' @param backfit_maxiter The maximum number of iterations used when backfitting the rank K fit
#' @param backfit_tol The convergence tolerance used for the backfit
#' @param optim_maxiter The maximum number of iterations used when optimizing each factor
#' @param optim_tol The convergence tolerance used for optimizing each factor
#'
#' @return A symEBcovMF object
#' @export
#'
sym_ebcovmf_backfit <- function(S, sym_ebcovmf_obj, ebnm_fn, backfit_maxiter = 100, backfit_tol = 10^(-8), optim_maxiter= 500, optim_tol = 10^(-8)){
  K <- length(sym_ebcovmf_obj$lambda)
  iter <- 1
  obj_diff <- Inf
  sym_ebcovmf_obj$backfit_vec_elbo_full <- NULL
  sym_ebcovmf_obj$backfit_iter_elbo_vec <- NULL

  # refit lambda
  # sym_ebcovmf_obj <- refit_lambda(S, sym_ebcovmf_obj, maxiter = 25)

  while((iter <= backfit_maxiter) && (obj_diff > backfit_tol)){
    # print(iter)
    obj_old <- sym_ebcovmf_obj$elbo
    # loop through each factor
    for (k in 1:K){
      # print(k)
      # compute residual matrix
      R <- S - tcrossprod(sym_ebcovmf_obj$L_pm[,-k, drop = FALSE] %*% diag(sqrt(sym_ebcovmf_obj$lambda[-k]), ncol = (K-1)))
      R2k <- compute_R2(S, sym_ebcovmf_obj$L_pm[,-k, drop = FALSE], sym_ebcovmf_obj$lambda[-k], (K-1)) #this is right but I have one instance where the values don't match what I expect

      # optimize factor
      factor_proposed <- optimize_factor(R, ebnm_fn, optim_maxiter, optim_tol, sym_ebcovmf_obj$L_pm[,k], sym_ebcovmf_obj$lambda[k], sym_ebcovmf_obj$fitted_gs[[k]], R2k, sym_ebcovmf_obj$n, sym_ebcovmf_obj$KL[-k])

      # update object
      # check if update leads to increase in objective function
      if ((factor_proposed$curr_elbo > sym_ebcovmf_obj$elbo) | (iter == 1)){
        sym_ebcovmf_obj$L_pm[,k] <- factor_proposed$v
        sym_ebcovmf_obj$KL[k] <- factor_proposed$rank_one_KL
        sym_ebcovmf_obj$lambda[k] <- factor_proposed$lambda_k
        sym_ebcovmf_obj$resid_s2 <- factor_proposed$resid_s2
        sym_ebcovmf_obj$fitted_gs[[k]] <- factor_proposed$fitted_g_k
        sym_ebcovmf_obj$elbo <- factor_proposed$curr_elbo
        sym_ebcovmf_obj$backfit_vec_elbo_full <- c(sym_ebcovmf_obj$backfit_vec_elbo_full, factor_proposed$vec_elbo_full)
      } else {
        obj_diff <- sym_ebcovmf_obj$elbo - factor_proposed$curr_elbo
        print(paste('update to factor', k, 'decreased the elbo by', abs(obj_diff)))
      }

      #print(sym_ebcovmf_obj$elbo)
      # sym_ebcovmf_obj <- refit_lambda(S, sym_ebcovmf_obj) # add refitting step?
      #print(sym_ebcovmf_obj$elbo)
    }
    sym_ebcovmf_obj$backfit_iter_elbo_vec <- c(sym_ebcovmf_obj$backfit_iter_elbo_vec, sym_ebcovmf_obj$elbo)

    iter <- iter + 1
    obj_diff <- abs(sym_ebcovmf_obj$elbo - obj_old)
  }
  # nullcheck
  sym_ebcovmf_obj <- nullcheck_factors(sym_ebcovmf_obj)
  return(sym_ebcovmf_obj)
}
