#' Optimize one factor
#'
#' @param R An n-by-n residual matrix
#' @param ebnm_fn An EBNM solver
#' @param maxiter The maximum number of iterations for the optimization
#' @param tol The convergence tolerance for the optimization
#' @param v_init An initial estimate for the factor
#' @param lambda_k An initial estimate for the factor weight, lambda
#' @param g_k_init An initial estimate for the prior; can be set to NULL
#' @param R2k The R2 value computed using the other factors
#' @param n The number of samples
#' @param KL The KL values for the other factors
#'
#' @return a proposed update for the factor estimate
#' @export
#'

optimize_factor <- function(R, ebnm_fn, maxiter, tol, v_init, lambda_k, g_k_init, R2k, n, KL){
  R2 <- R2k + lambda_k^2 - 2*lambda_k*as.numeric(t(v_init)%*% R %*%matrix(v_init, ncol = 1))
  resid_s2 <- estimate_resid_s2(n = n, R2 = R2)
  rank_one_KL <- 0
  curr_elbo <- -Inf
  obj_diff <- Inf
  fitted_g_k <- g_k_init
  iter <- 1
  vec_elbo_full <- NULL
  v <- v_init

  while((iter <= maxiter) && (obj_diff > tol)){
    # update l; power iteration step
    v.old <- v
    x <- R %*% v
    e <- ebnm_fn(x = x, s = sqrt(resid_s2), g_init = fitted_g_k)
    scaling_factor <- sqrt(sum(e$posterior$mean^2) + sum(e$posterior$sd^2))
    if (scaling_factor == 0){ # check if scaling factor is zero
      scaling_factor <- Inf
      v <- e$posterior$mean/scaling_factor
      lambda_k <- 0
      R2 <- R2k
      fitted_g_k <- e$fitted_g
      rank_one_KL <- as.numeric(e$log_likelihood) +
        - normal_means_loglik(x, sqrt(resid_s2), e$posterior$mean, e$posterior$mean^2 + e$posterior$sd^2)
      resid_s2 <- estimate_resid_s2(n = n, R2 = R2)
      curr_elbo <- compute_elbo(resid_s2 = resid_s2,
                                n = n,
                                KL = c(KL, rank_one_KL),
                                R2 = R2)
      print('Warning: scaling factor is zero')
      break
    }
    v <- e$posterior$mean/scaling_factor

    # update lambda and R2
    lambda_k.old <- lambda_k
    lambda_k <- max(as.numeric(t(v) %*% R %*% v), 0)
    R2 <- R2k - lambda_k^2

    #store estimate for g
    fitted_g_k.old <- fitted_g_k
    fitted_g_k <- e$fitted_g

    # store KL
    rank_one_KL.old <- rank_one_KL
    rank_one_KL <- as.numeric(e$log_likelihood) +
      - normal_means_loglik(x, sqrt(resid_s2), e$posterior$mean, e$posterior$mean^2 + e$posterior$sd^2)

    # update resid_s2
    resid_s2.old <- resid_s2
    resid_s2 <- estimate_resid_s2(n = n, R2 = R2) # this goes negative?????

    # check convergence - maybe change to rank-one obj function
    curr_elbo.old <- curr_elbo
    curr_elbo <- compute_elbo(resid_s2 = resid_s2,
                              n = n,
                              KL = c(KL, rank_one_KL),
                              R2 = R2)
    if (iter > 1){
      obj_diff <- curr_elbo - curr_elbo.old
    }
    if (obj_diff < 0){ # check if convergence_val < 0
      v <- v.old
      resid_s2 <- resid_s2.old
      rank_one_KL <- rank_one_KL.old
      lambda_k <- lambda_k.old
      curr_elbo <- curr_elbo.old
      fitted_g_k <- fitted_g_k.old
      print(paste('elbo decreased by', abs(obj_diff)))
      break
    }
    vec_elbo_full <- c(vec_elbo_full, curr_elbo)
    iter <- iter + 1
  }
  return(list(v = v, lambda_k = lambda_k, resid_s2 = resid_s2, curr_elbo = curr_elbo, vec_elbo_full = vec_elbo_full, fitted_g_k = fitted_g_k, rank_one_KL = rank_one_KL))
}
