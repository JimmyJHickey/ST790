library(jhickeyST790)

set.seed(12345)

n = 1e2
x = seq(0, 5, length.out = n)
y = sin(pi*x) + x + 0.5*rnorm(n)


# theta = 0, lambda = 100

theta = rep(0, n)
lambda = 10000

Dkn = myGetDkn(5, n)

L = hp_lipschitz_constant(lambda = lambda, Dkn = Dkn)

step_size = 1/L

hp_objective(y = y, theta = theta, Dkn = Dkn, lambda = lambda)
hp_objective_gradient(y = y, theta = theta, Dkn = Dkn, lambda = lambda)

x2 = gradient_step(gradf = hp_objective_gradient(y = y, theta = theta, Dkn = Dkn, lambda = lambda),
                   x = x, t = step_size)

x3 = gradient_step(gradf = hp_objective_gradient(y = y, theta = x2, Dkn = Dkn, lambda = lambda),
                   x = x2, t = step_size)

x4 = gradient_step(gradf = hp_objective_gradient(y = y, theta = x3, Dkn = Dkn, lambda = lambda),
                   x = x3, t = step_size)

hp_objective(y = y, theta = theta, Dkn = Dkn, lambda = lambda)
hp_objective(y = y, theta = x2, Dkn = Dkn, lambda = lambda)
hp_objective(y = y, theta = x3, Dkn = Dkn, lambda = lambda)
hp_objective(y = y, theta = x4, Dkn = Dkn, lambda = lambda)



# theta = y, lambda = 100

theta = y
lambda = 100
