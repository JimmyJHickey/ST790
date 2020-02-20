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



#' MM algorithm for smooth LAD regression
#'
#' @param y response
#' @param X design matrix
#' @param beta Initial regression coefficient vector
#' @param epsilon smoothing parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
smLAD <- function(y, X, beta, epsilon=0.25, max_iter=1e2, tol=1e-3) {

  # create vectors
  objective_history = c()
  relative_objective_history = c()
  relative_iterate_history = c()

  # initialize variables
  current_beta = beta
  objective_history[1] = fx_lad(y, X, current_beta, epsilon)
  relative_objective_history[1] = 0
  relative_iterate_history[1] = 0

  # perform mm update until either
  #   we have changed less than the tolerance
  #   we have done the maximum number of iterations
  for (i in 2:max_iter)
  {

    new_beta = mm_update(y, X, current_beta, epsilon = epsilon)

    objective_history[i] = fx_lad(y, X, new_beta, epsilon)

    relative_objective_history[i] = abs((objective_history[i] - objective_history[i-1]))/(1 + abs(objective_history[i]))
    relative_iterate_history[i] = norm(new_beta - current_beta, '2') / (1 + norm(new_beta, '2'))

    current_beta = new_beta

    # break if change less than tolerated amount
    if (relative_objective_history[i] <= tol)
      break

  } # end for

  return_list = list(
    "final_iterate" = current_beta,
    "objective_history" = objective_history,
    "relative_objective_history" = relative_objective_history,
    "relative_iterate_history" = relative_iterate_history
  )

  return(return_list)
}
