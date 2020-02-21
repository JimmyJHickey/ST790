#' Compute Newton Step (Naive) for logistic ridge regression
#'
#' @param y response
#' @param X Design matrix
#' @param beta Current regression vector estimate
#' @param g Gradient vector
#' @param lambda Regularization parameter
#' @param epsilon smoothing parameter
newton_step_naive <- function(y, X, beta, g, lambda, epsilon=0.25) {

  library(base)

  wBeta = diag( c( epsilon / ((y - X %*% beta )^2 + epsilon)^(3/2) ) )
  lhs = diag(lambda, ncol(X)) + t(X) %*% wBeta %*% X

  chol_out = chol(lhs)

  Ux = forwardsolve( t(chol_out), g )

  delta_beta_nt = -backsolve( chol_out, Ux )

  return(delta_beta_nt)
}



#' Compute Newton Step (Sherman-Morrison-Woodbury) for logistic ridge regression
#'
#'
#' @param y response
#' @param X Design matrix
#' @param beta Current regression vector estimate
#' @param g Gradient vector
#' @param lambda Regularization parameter
#' @param epsilon smoothing parameter
newton_step_smw <- function(y, X, beta, g, lambda, epsilon=0.25) {
  # (A + U C V)^-1 = A^-1 - A^-1 U ( C^-1 V A^-1 U)^-1 V A^-1

  # U, V^T
  X_wBeta_sqrt = t(X) %*% diag( c( sqrt(epsilon / ((y - X %*% beta )^2 + epsilon)^(3/2) ) ))

  # A
  lambda_diag = diag(1/lambda, ncol(X))

  # C
  id_p = diag(1, nrow(X))

  inv = lambda_diag - lambda_diag %*% X_wBeta_sqrt %*% solve(id_p + t(X_wBeta_sqrt) %*% lambda_diag %*% X_wBeta_sqrt) %*% t(X_wBeta_sqrt) %*% lambda_diag


  delta_beta_nt = -inv %*% g

  return(delta_beta_nt)
}
