#' Refit lambda
#'
#' @param S An n-by-n Gram/covariance matrix
#' @param sym_ebcovmf_obj A symEBcovMF object
#' @param maxiter The maximum number of iterations for optimization
#' @param tol The convergence tolerance parameter for optimization
#'
#' @return A symEBcovMF object
#' @export
#'
refit_lambda <- function(S, sym_ebcovmf_obj, maxiter = 100, tol = 10^(-6)){
  K <- length(sym_ebcovmf_obj$lambda)
  if (K <= 1){
    print('Cannot refit lambda')
  } else {
    sym_ebcovmf_obj.old <- sym_ebcovmf_obj
    iter <- 1
    obj_diff <- Inf
    curr_elbo <- -Inf
    while((iter <= maxiter) && (obj_diff > tol)){
      # print(iter)
      # Update lambdas
      lambda.old <- sym_ebcovmf_obj$lambda
      for (k in 1:K){
        Rk <- S - tcrossprod(sym_ebcovmf_obj$L_pm[,-k] %*% diag(sqrt(sym_ebcovmf_obj$lambda[-k]), ncol = (K-1)))
        sym_ebcovmf_obj$lambda[k] <- max(t(sym_ebcovmf_obj$L_pm[,k, drop = FALSE]) %*% Rk %*% sym_ebcovmf_obj$L_pm[,k, drop = FALSE], 0)
      }

      # Update resid_s2
      resid_s2.old <- sym_ebcovmf_obj$resid_s2
      sym_ebcovmf_obj$resid_s2 <- estimate_resid_s2(S = S,
                                    L = sym_ebcovmf_obj$L_pm,
                                    lambda = sym_ebcovmf_obj$lambda,
                                    n = sym_ebcovmf_obj$n,
                                    K = K)

      # Update elbo
      curr_elbo.old <- curr_elbo
      curr_elbo <- compute_elbo(S = S,
                                L = sym_ebcovmf_obj$L_pm,
                                lambda = sym_ebcovmf_obj$lambda,
                                resid_s2 = sym_ebcovmf_obj$resid_s2,
                                n = sym_ebcovmf_obj$n,
                                K = K,
                                KL = sym_ebcovmf_obj$KL)

      # Check convergence
      if (iter > 1){
        obj_diff <- curr_elbo - curr_elbo.old
      }
      if (obj_diff < 0){
        sym_ebcovmf_obj$lambda <- lambda.old
        sym_ebcovmf_obj$resid_s2 <- resid_s2.old
        curr_elbo <- curr_elbo.old
        print(paste('elbo decreased by', abs(obj_diff)))
        break
      }
      iter <- iter + 1
    }
    # nullcheck
    if (any(sym_ebcovmf_obj$lambda == 0)){
      idx <- which(sym_ebcovmf_obj$lambda != 0)
      sym_ebcovmf_obj$lambda <- sym_ebcovmf_obj$lambda[idx]
      sym_ebcovmf_obj$L_pm <- sym_ebcovmf_obj$L_pm[,idx]
      sym_ebcovmf_obj$KL <- sym_ebcovmf_obj$KL[idx]
      sym_ebcovmf_obj$resid_s2 <- estimate_resid_s2(S = S,
                                                    L = sym_ebcovmf_obj$L_pm,
                                                    lambda = sym_ebcovmf_obj$lambda,
                                                    n = sym_ebcovmf_obj$n,
                                                    K = K)
      curr_elbo <- compute_elbo(S = S,
                                L = sym_ebcovmf_obj$L_pm,
                                lambda = sym_ebcovmf_obj$lambda,
                                resid_s2 = sym_ebcovmf_obj$resid_s2,
                                n = sym_ebcovmf_obj$n,
                                K = K,
                                KL = sym_ebcovmf_obj$KL)
    }
    # check objective function
    if ((curr_elbo - sym_ebcovmf_obj.old$elbo) < 0){
      sym_ebcovmf_obj <- sym_ebcovmf_obj.old
    } else {
      sym_ebcovmf_obj$elbo <- curr_elbo
      sym_ebcovmf_obj$vec_elbo_K[K] <- curr_elbo
    }
  }
  return(sym_ebcovmf_obj)
}
