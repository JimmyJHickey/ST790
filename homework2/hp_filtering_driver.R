library(jhickeyST790)

set.seed(12345)

n = 1e2
x = seq(0, 5, length.out = n)
y = sin(pi*x) + x + 0.5*rnorm(n)
k = 5

# theta = 0, lambda = 100

theta = rep(0, n)
lambda = 10000

Dkn = myGetDkn(k, n)

L = hp_lipschitz_constant(lambda = lambda, Dkn = Dkn)
lipschitz_step_size = 1/L

hp_objective_wrapper = function(x)
  fx_hp(y = y, theta = x, Dkn = Dkn, lambda = lambda)

hp_gradient_wrapper = function(x)
  gradf_hp(y = y, theta = x, Dkn = Dkn, lambda = lambda)


gradient_descent_fixed( fx = hp_objective_wrapper,
                        gradf = hp_gradient_wrapper,
                        x0 = x,
                        t = lipschitz_step_size,
                        max_iter = 1e2,
                        tol = 1e-3
                        )
# testing


gradient_value = hp_gradient_wrapper(x)
step = gradient_step(gradient_value, x, lipschitz_step_size)
q = hp_objective_wrapper(step)

# hp_objective(y = y, theta = theta, Dkn = Dkn, lambda = lambda)
# hp_objective_gradient(y = y, theta = theta, Dkn = Dkn, lambda = lambda)
#
# x2 = gradient_step(gradf = hp_objective_gradient(y = y, theta = theta, Dkn = Dkn, lambda = lambda),
#                    x = x, t = step_size)
#
# x3 = gradient_step(gradf = hp_objective_gradient(y = y, theta = x2, Dkn = Dkn, lambda = lambda),
#                    x = x2, t = step_size)
#
# x4 = gradient_step(gradf = hp_objective_gradient(y = y, theta = x3, Dkn = Dkn, lambda = lambda),
#                    x = x3, t = step_size)
#
# hp_objective(y = y, theta = theta, Dkn = Dkn, lambda = lambda)
# hp_objective(y = y, theta = x2, Dkn = Dkn, lambda = lambda)
# hp_objective(y = y, theta = x3, Dkn = Dkn, lambda = lambda)
# hp_objective(y = y, theta = x4, Dkn = Dkn, lambda = lambda)



# theta = y, lambda = 100

theta = y
lambda = 100
