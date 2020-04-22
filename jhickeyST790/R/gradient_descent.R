#' Gradient Step
#'
#' @param gradf handle to function that returns gradient of objective function
#' @param x current parameter estimate
#' @param t step-size
#' @export
gradient_step <- function(gradf, x, t) {
  return(x - t * gradf)
}



#' Backtracking
#'
#' @param fx handle to function that returns objective function values
#' @param x current parameter estimate
#' @param t current step-size
#' @param df the value of the gradient of objective function evaluated at the current x
#' @param alpha the backtracking parameter
#' @param beta the decrementing multiplier
#' @export
backtrack <- function(fx, x, t, df, alpha=0.5, beta=0.9) {
  while( (fx(x - t * df ) >= fx(x) - t * alpha * norm(df, "2")^2)  & (t > 1e-10))
      t <- t * beta

  return(t)
}


#' Gradient Descent (Fixed Step-Size)
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function
#' @param x0 initial parameter estimate
#' @param t step-size
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @return A list with letters and numbers.
#' \itemize{
#'   \item{final_iterate}{The final iterate}
#'   \item{objective_history}{A vector of objective function value from each iteration.}
#'   \item{gradient_history}{A vector of the 2-norm of the gradient value from each iteration.}
#'   \item{relative_objective_history}{A vector of relative change in object value between iterations.}
#'   \item{relative_iterate_history}{A vector of relative change in iterate value between iterations.}
#' }
#' @export
gradient_descent_fixed <- function(fx, gradf, x0, t, max_iter=1e2, tol=1e-3) {

  # Check that step size is positive
  if (t <= 0)
    stop("Step size must be positive")

  # create vectors
  objective_history = c()
  gradient_history = c()
  relative_objective_history = c()
  relative_iterate_history = c()

  # initialize variables
  current_iterate = x0
  gradient_history[1] = 0
  objective_history[1] = fx(current_iterate)
  relative_objective_history[1] = 0
  relative_iterate_history[1] = 0

  # perform gradient descent until either
  #   we have changed less than the tolerance
  #   we have done the maximum number of iterations
  for (i in 2:max_iter)
  {

    # Calculate gradient for current x
    gradient_value = gradf(current_iterate)
    gradient_history[i] = norm(gradient_value, '2')

    # Gradient step to get new objective iterate value
    new_iterate = gradient_step(gradient_value, current_iterate, t)

    objective_history[i] = fx(new_iterate)

    relative_objective_history[i] = abs((objective_history[i] - objective_history[i-1]))/(1 + abs(objective_history[i]))
    relative_iterate_history[i] = norm(new_iterate - current_iterate, '2') / (1 + norm(new_iterate, '2'))

    current_iterate = new_iterate

    # break if change less than tolerated amount
    if (relative_objective_history[i] <= tol)
      break

  } # end for

  return_list = list(
    "final_iterate" = current_iterate,
    "objective_history" = objective_history,
    "gradient_history" = gradient_history,
    "relative_objective_history" = relative_objective_history,
    "relative_iterate_history" = relative_iterate_history
  )

  return(return_list)
}


#' Gradient Descent (Backtracking Step-Size)
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function
#' @param x0 initial parameter estimate
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent_backtrack <- function(fx, gradf, x0, max_iter=1e2, tol=1e-3)
{
  # create vectors
  objective_history = c()
  gradient_history = c()
  relative_objective_history = c()
  relative_iterate_history = c()

  # guess at step size
  default_step_size = 1

  # initialize variables
  current_iterate = x0
  gradient_history[1] = 0
  objective_history[1] = fx(current_iterate)
  relative_objective_history[1] = 0
  relative_iterate_history[1] = 0

  # perform gradient descent until either
  #   we have changed less than the tolerance
  #   we have done the maximum number of iterations
  for (i in 2:max_iter)
  {

    # Calculate gradient for current x
    gradient_value = gradf(current_iterate)
    gradient_history[i] = norm(gradient_value, '2')

    # Calculate step size using backtracking
    step_size = backtrack(fx = fx, x = current_iterate, t = default_step_size, df = gradient_value, alpha = 0.5, beta = 0.9)

    # Gradient step to get new objective iterate value
    new_iterate = gradient_step(gradient_value, current_iterate, step_size)

    objective_history[i] = fx(new_iterate)

    relative_objective_history[i] = abs((objective_history[i] - objective_history[i-1]))/(1 + abs(objective_history[i]))
    relative_iterate_history[i] = norm(new_iterate - current_iterate, '2') / (1 + norm(new_iterate, '2'))

    current_iterate = new_iterate

    # break if change less than tolerated amount
    if (relative_objective_history[i] <= tol)
      break

  } # end for

  return_list = list(
    "final_iterate" = current_iterate,
    "objective_history" = objective_history,
    "gradient_history" = gradient_history,
    "relative_objective_history" = relative_objective_history,
    "relative_iterate_history" = relative_iterate_history
  )

  return(return_list)
}


#' Gradient Descent Wrapper
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function
#' @param x0 initial parameter estimate
#' @param t step-size
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent <- function(fx, gradf, x0, t=NULL, max_iter=1e2, tol=1e-3) {

  if(is.null(t))
    return(gradient_descent_backtrack(fx = fx, gradf = gradf, x0 = x0, max_iter = max_iter, tol = tol))
  else
    return(gradient_descent_fixed(fx = fx, gradf = gradf, x0 = x0, t = t, max_iter = max_iter, tol = tol))
}
