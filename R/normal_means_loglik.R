#' Compute the normal means log likelihood
#'
#' @param x A data vector
#' @param s A vector of standard errors
#' @param Et A vector of posterior means
#' @param Et2 A vector of posterior second moments
#'
#' @return The log likelihood
#' 
#' @references This function is copied from flashier:::normal.means.loglik
#' (https://github.com/willwerscheid/flashier/blob/master/R/objective.R)
#' @export
#' 
normal_means_loglik <- function(x, s, Et, Et2) {
  idx <- is.finite(s) & s > 0
  x <- x[idx]
  s <- s[idx]
  Et <- Et[idx]
  Et2 <- Et2[idx]
  
  return(-0.5 * sum(log(2 * pi * s^2) + (1 / s^2) * (Et2 - 2 * x * Et + x^2)))
}