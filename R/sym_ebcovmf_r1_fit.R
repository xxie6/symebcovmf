#' Fit rank one symEBcovMF model
#'
#' @param S An n-by-n Gram/covariance matrix
#' @param sym_ebcovmf_obj A symEBcovMF object
#' @param ebnm_fn An EBNM solver
#' @param maxiter The maximum number of iterations used during optimization
#' @param tol The convergence tolerance parameter used during optimization
#' @param v_init An initial estimate for the loadings vector estimate. Default is NULL
#' @param sign_constraint The sign constraint (if any) for the loadings estimate
#'
#' @return A symEBcovMF object
#' @export
#'
sym_ebcovmf_r1_fit <- function(S, sym_ebcovmf_obj, ebnm_fn, maxiter, tol, v_init = NULL, sign_constraint = NULL){
  if (is.null(sym_ebcovmf_obj$L_pm) == FALSE){
    K <- length(sym_ebcovmf_obj$lambda) + 1
    R <- S - tcrossprod(sym_ebcovmf_obj$L_pm %*% diag(sqrt(sym_ebcovmf_obj$lambda), ncol = (K-1)))
    R2k <- compute_R2(S, sym_ebcovmf_obj$L_pm, sym_ebcovmf_obj$lambda, (K-1))
  } else {
    K <- 1
    R <- S
    R2k <- sum(S^2)
    sym_ebcovmf_obj$failyn <- FALSE
  }
  sym_ebcovmf_obj.old <- sym_ebcovmf_obj

  # initialize estimate for l
  if (is.null(v_init) == TRUE){
    # checks non-negativity constraint
    if (is.null(sign_constraint) == TRUE){
      nonnegative_constraint <- FALSE
    } else {
      nonnegative_constraint <- TRUE
    }
    sym_ebcovmf_v_init <- sym_ebcovmf_r1_init(R, nonnegative = nonnegative_constraint)
    v <- sym_ebcovmf_v_init$v
    lambda_k <- sym_ebcovmf_v_init$lambda_k
  } else {
    v <- v_init
    v <- v/sqrt(sum(v^2))
    lambda_k <- drop(t(v) %*% R %*% v)
  }

  # check if initial value of lambda is positive. if not, stop function.
  if (lambda_k <= 0){
    warning('Initial value of scaling factor is zero or nonnegative.')
    return(sym_ebcovmf_obj)
  }

  # initialize other values
  R2 <- R2k - lambda_k^2
  resid_s2 <- estimate_resid_s2(n = sym_ebcovmf_obj$n, R2 = R2)
  rank_one_KL <- 0
  curr_elbo <- -Inf
  obj_diff <- Inf
  fitted_g_k <- NULL
  iter <- 1

  sym_ebcovmf_obj$vec_elbo_full <- c(sym_ebcovmf_obj$vec_elbo_full, K)
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
      resid_s2 <- estimate_resid_s2(n = sym_ebcovmf_obj$n, R2 = R2)
      curr_elbo <- compute_elbo(resid_s2 = resid_s2,
                                n = sym_ebcovmf_obj$n,
                                KL = c(sym_ebcovmf_obj$KL, rank_one_KL),
                                R2 = R2)
      warning('scaling factor is zero')
      break
    }
    v <- e$posterior$mean/scaling_factor

    #store estimate for g
    fitted_g_k.old <- fitted_g_k
    fitted_g_k <- e$fitted_g

    # store KL
    rank_one_KL.old <- rank_one_KL
    rank_one_KL <- as.numeric(e$log_likelihood) +
      - normal_means_loglik(x, sqrt(resid_s2), e$posterior$mean, e$posterior$mean^2 + e$posterior$sd^2)

    R2 <- R2k + lambda_k^2 - 2*lambda_k*as.numeric(t(v)%*% R %*%v)
    ebnm_updated_elbo <- compute_elbo(resid_s2 = resid_s2,
                                      n = sym_ebcovmf_obj$n,
                                      KL = c(sym_ebcovmf_obj$KL, rank_one_KL),
                                      R2 = R2)

    # check 1: check if the EBNM problem leads to a large objective function decrease
    ebnm_updated_elbo_diff <- ebnm_updated_elbo - curr_elbo
    if (ebnm_updated_elbo_diff < 0){
      print(paste('Update to g and q of factor', K, 'decreased elbo by', abs(ebnm_updated_elbo_diff)))
      if (abs(ebnm_updated_elbo_diff) > 10^(-4)){
        sym_ebcovmf_obj$failyn <- TRUE
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
      print(paste('Update to g and q of factor', K, 'made v\'Rv negative'))
      sym_ebcovmf_obj$failyn <- TRUE
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
    resid_s2 <- estimate_resid_s2(n = sym_ebcovmf_obj$n, R2 = R2)

    # check if the residual variance estimate is negative
    # this shouldn't happen but I kept the check in to confirm empirically that this doesn't happen
    if (resid_s2 < 0){
      sym_ebcovmf_obj$failyn <- TRUE
      v <- v.old
      resid_s2 <- resid_s2.old
      rank_one_KL <- rank_one_KL.old
      lambda_k <- lambda_k.old
      curr_elbo <- curr_elbo.old
      fitted_g_k <- fitted_g_k.old
      print(paste('Update to factor', K, 'made the residual variance estimate negative'))
      break
    }

    # check convergence
    curr_elbo.old <- curr_elbo
    curr_elbo <- compute_elbo(resid_s2 = resid_s2,
                              n = sym_ebcovmf_obj$n,
                              KL = c(sym_ebcovmf_obj$KL, rank_one_KL),
                              R2 = R2)
    if (iter > 1){
      obj_diff <- curr_elbo - curr_elbo.old
    }
    if (obj_diff < 0){ # check if convergence_val < 0
      if (abs(obj_diff) > 10^(-4)){ # check if the objective function decrease is very large.
        sym_ebcovmf_obj$failyn <- TRUE
      }
      v <- v.old
      resid_s2 <- resid_s2.old
      rank_one_KL <- rank_one_KL.old
      lambda_k <- lambda_k.old
      curr_elbo <- curr_elbo.old
      fitted_g_k <- fitted_g_k.old
      print(paste('Update to factor', K, 'decreased elbo by', abs(obj_diff)))
      break
    }
    sym_ebcovmf_obj$vec_elbo_full <- c(sym_ebcovmf_obj$vec_elbo_full, curr_elbo)
    iter <- iter + 1
  }

  # nullcheck
  if((lambda_k == 0) | (sqrt(sum(v^2)) < 10^(-8))){ # if either of these hold, then that means the factor doesn't help the fit
    sym_ebcovmf_obj <- sym_ebcovmf_obj.old
  } else {
    sym_ebcovmf_obj$L_pm <- cbind(sym_ebcovmf_obj$L_pm, v)
    sym_ebcovmf_obj$KL[K] <- rank_one_KL
    sym_ebcovmf_obj$lambda[K] <- lambda_k
    sym_ebcovmf_obj$resid_s2 <- resid_s2
    sym_ebcovmf_obj$fitted_gs[[K]] <- fitted_g_k
    sym_ebcovmf_obj$elbo <- curr_elbo
    sym_ebcovmf_obj$vec_elbo_K <- c(sym_ebcovmf_obj$vec_elbo_K, curr_elbo)
  }
  return(sym_ebcovmf_obj)
}
