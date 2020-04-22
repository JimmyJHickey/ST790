#' Compute Lasso primal objective
#'
#' @param y response
#' @param b primal variable
#' @param D differencing matrix
#' @param lambda regularization parameter
#' @export
lasso_primal = function(y, b, D, lambda )
{
 return( 1/2 * norm(y - b, "2")^2 + lambda * norm( D %*% b, '1') )
}


#' Compute Lasso dual objective
#'
#' @param y response
#' @param v dual variable
#' @param D differencing matrix
#' @export
lasso_dual = function(y, v, D)
{
 return(1/2 * norm(y - t(D) %*% v, '2')^2 )
}


#' Compute Lasso dual objective original
#'
#' @param y response
#' @param v dual variable
#' @param D differencing matrix
#' @export
lasso_dual_original = function(y, v, D)
{
  return( 1/2 * norm(y, '2')^2 - 1/2 * norm(y - t(D) %*% v, '2')^2 )
}

#' Compute Lasso dual objective
#'
#' @param y response
#' @param v dual variable
#' @param D differencing matrix
#' @export
lasso_dual_gradient = function(y, v, D)
{
 return( D %*% t(D) %*% v - D %*% y )
}

#' Compute the primal variable b from the dual variable
#'
#' @param y response
#' @param v dual variable
#' @param D differencing matrix
#' @return the primal variable b
#' @export
dual_primal_map_b = function(y, v, D)
{
 return(y - t(D) %*% v)
}


#' Compute the primal variable theta from the dual variable
#'
#' @param y response
#' @param v dual variable
#' @param D differencing matrix
#' @return the primal variable b
#' @export
dual_primal_map_theta = function(y, v, D)
{
 return(D %*% dual_primal_map_b(y, v, D))
}

#' Compute duality gap
#'
#' @param y response
#' @param b primal variable
#' @param v dual variable
#' @param D differencing matrix
#' @param lambda regularization parameter
#' @export
duality_gap <- function(y, b, v, D, lambda) {
 return( lasso_primal(y, b, D, lambda) -lasso_dual_original(y, v, D) )
}


#' Compute KKT residual
#'
#' @param y response
#' @param b primal variable
#' @param theta primal variable
#' @param v dual variable
#' @param D differencing matrix
#' @param lambda regularization parameter
#' @export
kkt_residual <- function(y, b, theta, v, D, lambda) {

 tau = v

 # round theta so that values near 0 are 0
 theta = round(theta, 7)

 tau[theta != 0 ] = abs( abs(v[theta != 0]) - lambda )
 tau[theta == 0 ] = max( abs(v[theta == 0]) - lambda, 0)

 return(max(tau))
}


#' Compute Lasso proxmap
#'
#' @param v dual variable
#' @param lambda regularization parameter
#' @export
lasso_proxmap = function(v, lambda){
 return( ifelse(v > lambda, lambda,
                ifelse(v < - lambda, -lambda, v) ))
}

#' Proximal Gradient Step
#'
#' @param proxmap handle to a function that returns the proximal map
#' @param gradf gradient of the objective at x
#' @param x current parameter estimate
#' @param t step-size
#' @param lambda regularization parameter
#' @export
proximal_gradient_step <- function(proxmap, gradf, x, t, lambda) {
 return( proxmap( x - t * gradf, lambda) )
}


#' Solve trend-filtering by proximal gradient on dual problem
#'
#' @param y response
#' @param k order of differencing matrix
#' @param v Initial dual variable
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @return
#' \itemize{
#'   \item{final_iterate}{The final iterate}
#'   \item{objective_history}{A vector of objective function value from each iteration.}
#'   \item{relative_objective_history}{A vector of relative change in object value between iterations.}
#'   \item{relative_iterate_history}{A vector of relative change in iterate value between iterations.}
#'   \item{kkt_residual_history}{A vector of KKT residual value between iterations.}
#'   \item{duality_gap_history}{A vector of duality gap value between iterations.}
#' }
trend_filter_pg <- function(y, k, v, lambda=0, max_iter=1e2, tol=1e-3) {

 # Create differencing matrix
 n = length(y)
 D = as.matrix(myGetDkn(k, n))

 # Get Lipschitz step size
 t = 1 / norm(D, '2')^2

 b = dual_primal_map_b(y, v, D)
 theta = dual_primal_map_theta(y, v, D)

 # create vectors
 objective_history = c()
 relative_objective_history = c()
 relative_iterate_history = c()
 kkt_residual_history = c()
 duality_gap_history = c()


 # initialize variables
 current_iterate = v
 objective_history[1] = lasso_dual(y, current_iterate, D)
 relative_objective_history[1] = 0
 relative_iterate_history[1] = 0
 kkt_residual_history[1] = kkt_residual(y = y, b = b,
                                        theta = theta, v= current_iterate, D = D, lambda = lambda)
 duality_gap_history[1] = duality_gap(y, b, v, D, lambda)

 # perform gradient descent until either
 #   we have changed less than the tolerance
 #   we have done the maximum number of iterations
 for (i in 2:max_iter)
 {

   # Calculate gradient for current x
   gradient_value = lasso_dual_gradient(y, v, D)

   # Gradient step to get new objective iterate value
   new_iterate = proximal_gradient_step(proxmap = lasso_proxmap,
                                        gradf = gradient_value,
                                        x = current_iterate,
                                        t = t,
                                        lambda = lambda)

  b = dual_primal_map_b(y, new_iterate, D)
  theta = dual_primal_map_theta(y, new_iterate, D)

  objective_history[i] = lasso_dual(y, new_iterate, D)

  relative_objective_history[i] = abs((objective_history[i] - objective_history[i-1]))/(1 + abs(objective_history[i]))
  relative_iterate_history[i] = norm(new_iterate - current_iterate, '2') / (1 + norm(new_iterate, '2'))

  kkt_residual_history[i] = kkt_residual(y ,b, theta, new_iterate, D, lambda)
  duality_gap_history[i] = duality_gap(y, b, new_iterate, D, lambda)

  current_iterate = new_iterate

  # break if change less than tolerated amount
  if (duality_gap_history[i] <= tol)
   break

 } # end for

 return_list = list(
  "final_iterate" = current_iterate,
  "objective_history" = objective_history,
  "relative_objective_history" = relative_objective_history,
  "relative_iterate_history" = relative_iterate_history,
  "kkt_residual_history" = kkt_residual_history,
  "duality_gap_history" = duality_gap_history
 )

 return(return_list)
}


#' Coordinate descent update
#'
#' @param y response
#' @param v dual variable
#' @param lambda regularization parameter
#' @export
coordinate_descent = function(v, lambda, D, y)
{
  v_update = v
  n = length(v)
  bool_array = rep(TRUE, n)

  for (i in 1:n)
  {
    bool_array_i = bool_array
    bool_array_i[i] = FALSE

    # sum( d_ij * y_j, j=1:n)
    numerator = sum( D[i, ] * y)

    # numerator - \sum( v_k * \sum( d_kj* j_ij, j=1:n)  , k != i)
    numerator = numerator + sum( v[bool_array_i] * D[bool_array_i, ] * D[i,] )

    # denominator = sum( d_ij^2 , j: 1, n)
    denominator = sum( D[i, ]^2 )

    v_update[i] = numerator / denominator

    v_update[i] = lasso_proxmap(v_update[i], lambda)
  }

  return(v_update)
}


#' Solve trend-filtering by coordinate descent on dual problem
#'
#' @param y response
#' @param k order of differencing matrix
#' @param v Initial dual variable
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @return
#' \itemize{
#'   \item{final_iterate}{The final iterate}
#'   \item{objective_history}{A vector of objective function value from each iteration.}
#'   \item{relative_objective_history}{A vector of relative change in object value between iterations.}
#'   \item{relative_iterate_history}{A vector of relative change in iterate value between iterations.}
#'   \item{kkt_residual_history}{A vector of KKT residual value between iterations.}
#'   \item{duality_gap_history}{A vector of duality gap value between iterations.}
#' }
trend_filter_cd <- function(y, k, v, lambda=0, max_iter=1e2, tol=1e-3) {

  # Create differencing matrix
  n = length(y)
  D = as.matrix(myGetDkn(k, n))

  # Get Lipschitz step size
  t = 1 / norm(D, '2')^2

  b = dual_primal_map_b(y, v, D)
  theta = dual_primal_map_theta(y, v, D)

  # create vectors
  objective_history = c()
  relative_objective_history = c()
  relative_iterate_history = c()
  kkt_residual_history = c()
  duality_gap_history = c()


  # initialize variables
  current_iterate = v
  objective_history[1] = lasso_dual(y, current_iterate, D)
  relative_objective_history[1] = 0
  relative_iterate_history[1] = 0
  kkt_residual_history[1] = kkt_residual(y = y, b = b,
                                         theta = theta, v= current_iterate, D = D, lambda = lambda)
  duality_gap_history[1] = duality_gap(y, b, v, D, lambda)

  # perform gradient descent until either
  #   we have changed less than the tolerance
  #   we have done the maximum number of iterations
  for (i in 2:max_iter)
  {

    # Coordinate Descent step
    new_iterate = coordinate_descent(v = current_iterate,
                                     lambda = lambda,
                                     D = D,
                                     y = y)

    b = dual_primal_map_b(y, new_iterate, D)
    theta = dual_primal_map_theta(y, new_iterate, D)

    objective_history[i] = lasso_dual(y, new_iterate, D)

    relative_objective_history[i] = abs((objective_history[i] - objective_history[i-1]))/(1 + abs(objective_history[i]))
    relative_iterate_history[i] = norm(new_iterate - current_iterate, '2') / (1 + norm(new_iterate, '2'))

    kkt_residual_history[i] = kkt_residual(y ,b, theta, new_iterate, D, lambda)
    duality_gap_history[i] = duality_gap(y, b, new_iterate, D, lambda)

    current_iterate = new_iterate

    # break if change less than tolerated amount
    if (duality_gap_history[i] <= tol)
      break

  } # end for

  return_list = list(
    "final_iterate" = current_iterate,
    "objective_history" = objective_history,
    "relative_objective_history" = relative_objective_history,
    "relative_iterate_history" = relative_iterate_history,
    "kkt_residual_history" = kkt_residual_history,
    "duality_gap_history" = duality_gap_history
  )

  return(return_list)
}
