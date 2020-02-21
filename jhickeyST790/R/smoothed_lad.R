#' Objective Function for ridge LAD regression
#'
#' @param y response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param epsilon smoothing parameter
#' @param lambda regularization parameter
#' @export
fx_lad <- function(y, X, beta, epsilon=0.25,lambda=0) {
  residual = (y - X %*% beta)^2
  objective = sum( sqrt(residual + epsilon)) + lambda/2 * norm(beta, '2')^2

  return(objective)
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
  residual = (y - X %*% beta)^2
  sum_diag = diag( as.vector(1 / sqrt(residual + epsilon) ) )
  numerator = ( X %*% beta - y)

  gradient = t(X) %*% sum_diag %*% numerator + lambda * beta

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


#' Damped Newton's Method for Fitting Ridge LAD Regression
#'
#' @param y response
#' @param X Design matrix
#' @param beta Initial regression coefficient vector
#' @param epsilon smoothing parameter
#' @param lambda regularization parameter
#' @param naive Boolean variable; TRUE if using Cholesky on the Hessian
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
lad_newton <- function(y, X, beta, epsilon=0.25,lambda=0, naive=TRUE, max_iter=1e2, tol=1e-3) {

  newton = if(naive) newton_step_naive else newton_step_smw
  fx_lad_wrapper = function(in_beta) { fx_lad(y=y, X=X, beta=in_beta, epsilon = epsilon, lambda = lambda  ) }
  current_beta = beta
  t0 = 1

  for (i in 1:max_iter)
  {
    gradient = gradf_lad(y = y, X = X, beta = current_beta, epsilon = epsilon, lambda = lambda)

    newton_step = newton(y = y, X = X, beta = current_beta, g = gradient, epsilon = epsilon, lambda = lambda)
    lambda_condition = t(gradient) %*% newton_step

    if(lambda_condition^2/2 <= epsilon)
        break

    step_size = backtrack_descent(fx=fx_lad_wrapper, x=current_beta, t=t0, df=gradient, d=newton_step, alpha=0.5, beta=0.9)

    current_beta = current_beta + step_size * newton_step
  }

  return(current_beta)

}
