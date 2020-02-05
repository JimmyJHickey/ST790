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
  while( fx(x - t * df ) >= fx(x) - 0.5 * t * alpha * (t(df) %*% df) )
    t = t * beta

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

  # initialize variables
  current_x = x0

  gradient_values = vector()
  gradient_values[1] = 0
  objective = vector()
  objective[1] = fx(current_x)
  steps = vector()
  steps[1]=0

  # perform gradient descent until either
  #   we have changed less than the tolerance
  #   we have done the maximum number of iterations
  for (i in 2:max_iter)
  {

    # Calculate gradient for current x
    # gradient_values[i] = gradf(current_x)


    grad_value = gradf(current_x)


    # Gradient step to get new objective function value
    current_x = grad_step = gradient_step(grad_value, current_x, t)

    objective[i] = fx(grad_step)

    # iterate_change[i] = gradient_values[i] - gradient_values[i-1]
    # objective_change[i] = abs(objective[i] - objective[i-1])/objective[i-1]

    # break if change less than tolerated amount
    # if (objective_change[i] <= tol)
    #   break

  } # end for

  return_list = list(
    "final_iterate" = current_x,
    "objective_values" = objective,
    # "gradient_values" = gradient_values,
    # "objective_change" = objective_change,
    "iterate_change" = 1
  )

  return(return_list)
}



