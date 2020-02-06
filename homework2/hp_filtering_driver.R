library(jhickeyST790)

set.seed(12345)

n = 1e2
x = seq(0, 5, length.out = n)
y = sin(pi*x) + x + 0.5*rnorm(n)
k = 1

# theta = 0-5, lambda = 100

lambda = 100

Dkn = myGetDkn(k, n)

L = hp_lipschitz_constant(lambda = lambda, Dkn = Dkn)
lipschitz_step_size = 1/L

hp_objective_wrapper = function(x)
  fx_hp(y = y, theta = x, Dkn = Dkn, lambda = lambda)

hp_gradient_wrapper = function(x)
  gradf_hp(y = y, theta = x, Dkn = Dkn, lambda = lambda)


###
# theta = 0, lambda = 100
###
theta=rep(0, n)

# Without backtracking
zero_no_backtracking = gradient_descent( fx = hp_objective_wrapper,
                  gradf = hp_gradient_wrapper,
                  x0 = theta,
                  t = lipschitz_step_size,
                  max_iter = 1e3,
                  tol = 1e-3
                  )
plot(zero_no_backtracking$objective_history - zero_no_backtracking$objective_history[length(zero_no_backtracking$objective_history)])

# With backtracking
zero_w_backtracking = gradient_descent( fx = hp_objective_wrapper,
                  gradf = hp_gradient_wrapper,
                  x0 = theta,
                  max_iter = 1e3,
                  tol = 1e-3
)
plot(zero_w_backtracking$objective_history - zero_w_backtracking$objective_history[length(zero_w_backtracking$objective_history)])


###
# theta = y, lambda = 100
###

theta = y
# Without backtracking
y_no_backtracking = gradient_descent( fx = hp_objective_wrapper,
                                         gradf = hp_gradient_wrapper,
                                         x0 = theta,
                                         t = lipschitz_step_size,
                                         max_iter = 1e3,
                                         tol = 1e-3
)
plot(y_no_backtracking$objective_history - y_no_backtracking$objective_history[length(y_no_backtracking$objective_history)])

# With backtracking
y_w_backtracking = gradient_descent( fx = hp_objective_wrapper,
                                        gradf = hp_gradient_wrapper,
                                        x0 = theta,
                                        max_iter = 1e3,
                                        tol = 1e-3
)
plot(y_w_backtracking$objective_history - y_w_backtracking$objective_history[length(y_w_backtracking$objective_history)])


###
# time series data
###
library(tseries)
apple_data = get.hist.quote( instrument = "AAPL",
                             start = "2018-01-01",
                             end = "2020-01-01",
                             quote="Close")
