###
#
# Jimmy Hickey
# 2020-03-05
#
# Run l1 trend filtering code
#
###

# apple closing prices 2018-2020
library(tseries)
library(ggplot2)

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

  l1_fit = myTrendFilter(y, k, lambda)$x[1:n]

  aapl_df = data.frame(cbind(time_stock, l1_fit))

  ggplot(data = aapl_df, aes(x = times)) +
   geom_point(aes(y = Close), size = 0.1) +
   geom_line(aes(y = l1_fit)) +
   ggtitle(paste("Apple Stock Prices: ", expression(lambda), " =", lambda , ", k = ", k ))
 }

apple_plots(k = 3, lambda = 0.5)
apple_plots(k = 3, lambda = 10)
apple_plots(k = 3, lambda = 15)

apple_plots(k = 7, lambda = 0.5)
apple_plots(k = 7, lambda = 10)
apple_plots(k = 7, lambda = 15)
