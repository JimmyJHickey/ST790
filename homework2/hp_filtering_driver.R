library(jhickeyST790)
library(ggplot2)

set.seed(12345)

n = 1e2
x = seq(0, 5, length.out = n)
y = sin(pi*x) + x + 0.5*rnorm(n)
k = 1

# theta = 0-5, lambda = 100

lambda = 100

setup = function(y, k, n, lambda)
{
  Dkn = myGetDkn(k, n)
  L = hp_lipschitz_constant(lambda = lambda, Dkn = Dkn)
  lipschitz_step_size = 1/L

  hp_objective_wrapper = function(x)
    fx_hp(y = y, theta = x, Dkn = Dkn, lambda = lambda)

  hp_gradient_wrapper = function(x)
    gradf_hp(y = y, theta = x, Dkn = Dkn, lambda = lambda)

  output_list = list(
    "Dkn" = Dkn,
    "step_size" = lipschitz_step_size,
    "hp_fx" = hp_objective_wrapper,
    "hp_grad" = hp_gradient_wrapper
  )

  return(output_list)
}

orig_setup = setup(y=y, k=k, n=n, lambda = lambda)

Dkn = orig_setup$Dkn
lipschitz_step_size = orig_setup$step_size
hp_objective_wrapper = orig_setup$hp_fx
hp_gradient_wrapper = orig_setup$hp_grad


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

qplot(seq_along(zero_no_backtracking$objective_history),
      zero_no_backtracking$objective_history - zero_no_backtracking$objective_history[length(zero_no_backtracking$objective_history)]) +
  xlab("iteration") +
  ylab(expression(l(theta[m]) - l(theta[max]))) +
  ggtitle(expression(x[0] == 0), "fixed step size")

zero_df = data.frame(cbind(x, y, zero_no_backtracking$final_iterate))
ggplot(data = zero_df, aes(x = x)) + geom_point(aes(y = y)) + geom_line(aes(y = V3)) + ggtitle(expression(x[0] == 0), "fixed step size")


# With backtracking
zero_w_backtracking = gradient_descent( fx = hp_objective_wrapper,
                  gradf = hp_gradient_wrapper,
                  x0 = theta,
                  max_iter = 1e3,
                  tol = 1e-3
)

qplot(seq_along(zero_w_backtracking$objective_history),
  zero_w_backtracking$objective_history - zero_w_backtracking$objective_history[length(zero_w_backtracking$objective_history)]) +
    xlab("iteration") +
  ylab(expression(l(theta[m]) - l(theta[max]))) +
  ggtitle(expression(x[0] == 0), "backtracking")



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

qplot(seq_along(y_no_backtracking$objective_history),
      y_no_backtracking$objective_history - y_no_backtracking$objective_history[length(y_no_backtracking$objective_history)]) +
  xlab("iteration") +
  ylab(expression(l(theta[m]) - l(theta[max]))) +
  ggtitle(expression(x[0] == y), "fixed")

y_df = data.frame(cbind(x, y, zero_no_backtracking$final_iterate))

ggplot(data = y_df, aes(x = x)) + geom_point(aes(y = y)) + geom_line(aes(y = V3))+ ggtitle(expression(x[0] == y)) +   ggtitle(expression(x[0] == y), "fixed")




# With backtracking
y_w_backtracking = gradient_descent( fx = hp_objective_wrapper,
                                        gradf = hp_gradient_wrapper,
                                        x0 = theta,
                                        max_iter = 1e3,
                                        tol = 1e-3
                                    )



qplot(seq_along(y_w_backtracking$objective_history),
      y_w_backtracking$objective_history - y_w_backtracking$objective_history[length(y_w_backtracking$objective_history)]) +
  xlab("iteration") +
  ylab(expression(l(theta[m]) - l(theta[max]))) +
  ggtitle(expression(x[0] == y), "backtracking")




###
# time series data
###
library(tseries)
apple_data = get.hist.quote( instrument = "AAPL",
                             start = "2018-01-01",
                             end = "2020-01-01",
                             quote="Close")

times = time(apple_data)
time_stock = cbind(times, as.data.frame(apple_data))


apple_plots = function(k, lambda)
{
  n = nrow(time_stock)

  y=time_stock[,2]
  aapl_setup = setup(y=y, k=k, n=n, lambda = lambda)


  Dkn = aapl_setup$Dkn
  lipschitz_step_size = aapl_setup$step_size
  hp_objective_wrapper = aapl_setup$hp_fx
  hp_gradient_wrapper = aapl_setup$hp_grad

  aapl_no_backtracking = gradient_descent( fx = hp_objective_wrapper,
                                           gradf = hp_gradient_wrapper,
                                           x0 = y,
                                           t = step_size,
                                           max_iter = 1e3,
                                           tol = 1e-3
  )

  aapl_1_df = data.frame(cbind(time_stock, aapl_no_backtracking$final_iterate))
  ggplot(data = aapl_1_df, aes(x = times)) +
    geom_point(aes(y = Close), size = 0.1) +
    geom_line(aes(y = aapl_no_backtracking.final_iterate)) +
    ggtitle(paste("Apple Stock Prices: ", expression(lambda), " =", lambda , ", k = ", k ))
}


apple_plots(k = 3, lambda = 0.005)
apple_plots(k = 3, lambda = 0.5)
apple_plots(k = 3, lambda = 15)


apple_plots(k = 5, lambda = 0.005)
apple_plots(k = 5, lambda = 0.5)
apple_plots(k = 5, lambda = 15)


apple_plots(k = 2, lambda = 100)
apple_plots(k = 3, lambda = 100)

