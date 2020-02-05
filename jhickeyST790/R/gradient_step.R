#' Gradient Step
#'
#' @param gradf handle to function that returns gradient of objective function
#' @param x current parameter estimate
#' @param t step-size
#' @export
gradient_step <- function(gradf, x, t) {
  return(x - t * gradf)
}


#' Gradient Descent (Fixed Step-Size)
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function
#' @param x0 initial parameter estimate
#' @param t step-size
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent_fixed <- function(fx, gradf, x0, t, max_iter=1e2, tol=1e-3) {

  # Check that step size is positive
  if (t <= 0)
    stop("Step size must be positive")


  # initialize variables
  current_x = x0
  gradient_values[1] = 0


  # perform gradient descent until either
  #   we have changed less than the tolerance
  #   we have done the maximum number of iterations
  for (i in 2:max_iter)
  {


    gradient_values[i] = gradf(y = y, theta = current_x, Dkn = Dkn, lambda = lambda)
    iterate_change[i] = gradient_values[i] - gradient_values[i-1]

    objective[i] = gradient_step(gradf(fx, current_x), current_x, t)
    current_x = objective[i]
    objective_change[i] = abs(objective[i] - objective[i-1])/objective[i-1]


    # break if change less than tolerated amount
    if (objective_change[i] <= tol)
      break

  } # end for

  return_list = list(
    "final_iterate" = current_x,
    "objective_values" = objective,
    "gradient_values" = gradient_values,
    "objective_change" = objective_change,
    "iterate_change" = 1
  )

  return(return_list)
}
