#' Fit rank one symEBcovMF model
#'
#' @param S An n-by-n Gram/covariance matrix
#' @param sym_ebcovmf_obj A symEBcovMF object
#' @param ebnm_fn An EBNM solver
#' @param maxiter The maximum number of iterations used during optimization
#' @param tol The convergence tolerance parameter used during optimization
#'
#' @return A symEBcovMF object
#' @export
#'
sym_ebcovmf_r1_fit <- function(S, sym_ebcovmf_obj, ebnm_fn, maxiter, tol){
  if (is.null(sym_ebcovmf_obj$L_pm) == FALSE){
    K <- length(sym_ebcovmf_obj$lambda) + 1
    R <- S - tcrossprod(sym_ebcovmf_obj$L_pm %*% diag(sqrt(sym_ebcovmf_obj$lambda), ncol = (K-1)))
    R2k <- compute_R2(S, sym_ebcovmf_obj$L_pm, sym_ebcovmf_obj$lambda, (K-1))
  } else {
    K <- 1
    R <- S
    R2k <- sum(S^2)
  }
  sym_ebcovmf_obj.old <- sym_ebcovmf_obj

  # initialize estimate for l
  sym_ebcovmf_v_init <- sym_ebcovmf_r1_init(R)
  v <- sym_ebcovmf_v_init$v
  lambda_k <- sym_ebcovmf_v_init$lambda_k

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
    resid_s2 <- estimate_resid_s2(n = sym_ebcovmf_obj$n, R2 = R2)

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
      v <- v.old
      resid_s2 <- resid_s2.old
      rank_one_KL <- rank_one_KL.old
      lambda_k <- lambda_k.old
      curr_elbo <- curr_elbo.old
      fitted_g_k <- fitted_g_k.old
      print(paste('elbo decreased by', abs(obj_diff)))
      break
    }
    sym_ebcovmf_obj$vec_elbo_full <- c(sym_ebcovmf_obj$vec_elbo_full, curr_elbo)
    iter <- iter + 1
  }

  # nullcheck
  if((lambda_k == 0) | (sqrt(sum(v^2)) < 10^(-8))){
    #print('additional factor does not improve fit')
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
