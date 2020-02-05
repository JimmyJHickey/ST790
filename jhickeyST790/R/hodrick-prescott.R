#' Hodrick-Prescott filtering objective function
#'
#' @param y response
#' @param theta regression coefficient vector
#' @param Dkn sparse differencing matrix
#' @param lambda regularization parameter
#' @export
hp_objective <- function(y, theta, Dkn, lambda = 0)
{
  return(0.5 * sum((y-theta)^2) + (lambda/2) * sum((as.vector(Dkn %*% theta))^2) )
}


#' Hodrick-Prescott filtering objective gradient
#'
#' @param y response
#' @param theta regression coefficient vector
#' @param Dkn sparse differencing matrix
#' @param lambda regularization parameter
#' @export
hp_objective_gradient <- function(y, theta, Dkn, lambda = 0)
{
  return(-y + theta + lambda * as.vector(crossprod(Dkn) %*% theta) )
}


#' Get Lipschitz constant for Hodrick-Prescott filtering
#'
#' @param lambda regularization parameter
#' @param Dkn sparse differencing matrix
hp_lipschitz_constant <- function(lambda, Dkn)
{
  Dkn_operator_norm = max(svd(Dkn)$d)^2
  L = 1 + lambda *Dkn_operator_norm

  return(L)
}
