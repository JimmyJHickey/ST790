#' Objective Function for ridge LAD regression
#'
#' @param y response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param epsilon smoothing parameter
#' @param lambda regularization parameter
#' @export
fx_lad <- function(y, X, beta, epsilon=0.25,lambda=0) {
  objective = 0

  for(i in 1:nrow(X))
    objective = objective + sqrt( (y[i] - t(X[i,]) %*% beta)^2 + epsilon)  + (lambda / 2) * norm(beta,'2')^2

  return(as.numeric(objective))
}

#' Gradient for ridge LAD regression
#'
#' @param y response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param epsilon smoothing parameter
#' @param lambda regularization parameter
#' @export
gradf_lad <- function(y, X, beta, epsilon=0.25, lambda=0) {
  gradient = 0

  for(i in 1:nrow(X))
    gradient = gradient + ( (y[i] - t(X[i,]) %*% beta)^2 + epsilon)^(-1/2) * (y[i]- t(X[i,]) %*% beta) %*% -X[i,] + lambda * norm(beta,'2')

  return(gradient)
}

#' MM update for LAD regression
#'
#' @param y response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param epsilon smoothing parameter
#' @export
mm_update = function(y, X, beta, epsilon=0.25){

  W_diags = as.vector(0)

  for(i in 1:nrow(X))
    W_diags[i] =   1/( sqrt( (y[i] - t(X[i,]) %*% beta)^2 + epsilon) )

  W = diag( W_diags )

  update = solve( t(X) %*% W %*% X, t(X) %*% W %*% y )
  return(update)
}
