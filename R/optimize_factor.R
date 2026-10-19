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
  # Initialize values
  R2 <- R2k + lambda_k^2 - 2*lambda_k*as.numeric(t(v_init)%*% R %*%matrix(v_init, ncol = 1))
  resid_s2 <- estimate_resid_s2(n = n, R2 = R2)
  rank_one_KL <- 0
  curr_elbo <- -Inf
  obj_diff <- Inf
  fitted_g_k <- g_k_init
  iter <- 1
  vec_elbo_full <- NULL
  v <- v_init
  failyn <- FALSE

  while((iter <= maxiter) && (obj_diff > tol)){
    # update l; power iteration step
    v.old <- v
    x <- R %*% v
    if (iter == 1){
      s2 <- resid_s2*(as.numeric(t(v) %*% R %*% v)/lambda_k)
    } else {
      s2 <- resid_s2
    }
    e <- ebnm_fn(x = x,
                 s = sqrt(s2),
                 g_init = fitted_g_k)
    scaling_factor <- sqrt(sum(e$posterior$mean^2) + sum(e$posterior$sd^2))
    if (scaling_factor == 0){ # check if scaling factor is zero
      scaling_factor <- Inf
      v <- e$posterior$mean/scaling_factor
      lambda_k <- 0
      R2 <- R2k
      fitted_g_k <- e$fitted_g
      rank_one_KL <- as.numeric(e$log_likelihood) +
        - normal_means_loglik(x, sqrt(s2), e$posterior$mean, e$posterior$mean^2 + e$posterior$sd^2)
      resid_s2 <- estimate_resid_s2(n = n, R2 = R2)
      curr_elbo <- compute_elbo(resid_s2 = resid_s2,
                                n = n,
                                KL = c(KL, rank_one_KL),
                                R2 = R2)
      print('Warning: scaling factor is zero')
      break
    }
    v <- e$posterior$mean/scaling_factor

    #store estimate for g
    fitted_g_k.old <- fitted_g_k
    fitted_g_k <- e$fitted_g

    # store KL
    rank_one_KL.old <- rank_one_KL
    rank_one_KL <- as.numeric(e$log_likelihood) +
      - normal_means_loglik(x, sqrt(s2), e$posterior$mean, e$posterior$mean^2 + e$posterior$sd^2)

    R2 <- R2k + lambda_k^2 - 2*lambda_k*as.numeric(t(v)%*% R %*%v)
    ebnm_updated_elbo <- compute_elbo(resid_s2 = resid_s2,
                                      n = n,
                                      KL = c(KL, rank_one_KL),
                                      R2 = R2)

    # check 1: check if the EBNM problem leads to a large objective function decrease
    ebnm_updated_elbo_diff <- ebnm_updated_elbo - curr_elbo
    if (ebnm_updated_elbo_diff < 0){
      print(paste('Update to g and q decreased elbo by', abs(ebnm_updated_elbo_diff)))
      if (abs(ebnm_updated_elbo_diff) > 10^(-4)){
        failyn <- TRUE
        # need to set estimate to old values
        v <- v.old
        fitted_g_k <- fitted_g_k.old
        rank_one_KL <- rank_one_KL.old
        R2 <- R2k - lambda_k^2
        break #exit the for loop if this holds true
      }
    }
    # check 2: check if v'Rv is negative
    if (as.numeric(t(v) %*% R %*% v) < 0){
      print('Update to g and q made v\'Rv negative')
      failyn <- TRUE
      v <- v.old
      fitted_g_k <- fitted_g_k.old
      rank_one_KL <- rank_one_KL.old
      R2 <- R2k - lambda_k^2
      break #exit the for loop if this holds true
    }

    # update lambda and R2
    lambda_k.old <- lambda_k
    lambda_k <- max(as.numeric(t(v) %*% R %*% v), 0)
    R2 <- R2k - lambda_k^2

    # update resid_s2
    resid_s2.old <- resid_s2
    resid_s2 <- estimate_resid_s2(n = n, R2 = R2)

    # check if the residual variance estimate is negative
    # this shouldn't happen but I kept the check in to confirm empirically that this doesn't happen
    if (resid_s2 < 0){
      failyn <- TRUE
      v <- v.old
      resid_s2 <- resid_s2.old
      rank_one_KL <- rank_one_KL.old
      lambda_k <- lambda_k.old
      curr_elbo <- curr_elbo.old
      fitted_g_k <- fitted_g_k.old
      print(paste('Updated residual variance estimate is negative'))
      break
    }

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
      if (abs(obj_diff) > 10^(-4)){
        failyn <- TRUE
      }
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
  return(list(v = v, lambda_k = lambda_k, resid_s2 = resid_s2, curr_elbo = curr_elbo, vec_elbo_full = vec_elbo_full, fitted_g_k = fitted_g_k, rank_one_KL = rank_one_KL, failyn = failyn))
}
