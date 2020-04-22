library(ggplot2)

library(tseries)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

make_plots = function(y, k, lambda, max_iter=1e2, tol=1e-3)
{
  n = length(y)
  v0 = rep(1, n-k)
  D = as.matrix(myGetDkn(k, n))

  cd_out = trend_filter_cd(y=y,
                           k=k,
                           v=v0,
                           lambda = lambda,
                           max_iter = max_iter,
                           tol = tol)

  cd_dual_out = y - t(D) %*% cd_out$final_iterate

  aapl_1_df = data.frame(cbind(time_stock, cd_dual_out))
  cd_plot = ggplot(data = aapl_1_df, aes(x = times)) +
    geom_point(aes(y = Close), size = 0.1) +
    geom_line(aes(y = cd_dual_out)) +
    ggtitle(paste("CD Apple Stock Prices: ", expression(lambda), " =", lambda , ", k = ", k ))

  plot(cd_plot)

  iter_plot = ggplot() +
    geom_point(aes( y= cd_out$relative_iterate_history, x = 1:length(cd_out$duality_gap_history))) +
    ggtitle(paste("CD Relative Iterate History: ", expression(lambda), " =", lambda , ", k = ", k )) +
    ylab("Relative Iterate") +
    xlab("Iteration")

  obj_plot = ggplot() +
    geom_point(aes( y= cd_out$relative_objective_history, x = 1:length(cd_out$relative_objective_history))) +
    ggtitle(paste("CD Relative Objective History: ", expression(lambda), " =", lambda , ", k = ", k )) +
    ylab("Relative Objective") +
    xlab("Iteration")

  dg_plot = ggplot() +
    geom_point(aes( y= cd_out$duality_gap_history, x = 1:length(cd_out$duality_gap_history))) +
    ggtitle(paste("CD Duality Gap History: ", expression(lambda), " =", lambda , ", k = ", k )) +
    ylab("Duality Gap") +
    xlab("Iteration")

  kkt_plot = ggplot() +
    geom_point(aes( y= cd_out$kkt_residual_history, x = 1:length(cd_out$duality_gap_history))) +
    ggtitle(paste("CD KKT Residual History: ", expression(lambda), " =", lambda , ", k = ", k )) +
    ylab("KKT Residual") +
    xlab("Iteration")

  multiplot(iter_plot, obj_plot, dg_plot, kkt_plot, cols=2)

  pg_out = trend_filter_pg(y=y,
                           k=k,
                           v=v0,
                           lambda = lambda,
                           max_iter = max_iter,
                           tol = tol)

  pg_dual_out = y - t(D) %*% pg_out$final_iterate

  aapl_1_df = data.frame(cbind(time_stock, pg_dual_out))
  pg_plot = ggplot(data = aapl_1_df, aes(x = times)) +
    geom_point(aes(y = Close), size = 0.1) +
    geom_line(aes(y = cd_dual_out)) +
    ggtitle(paste("PG Apple Stock Prices: ", expression(lambda), " =", lambda , ", k = ", k ))

  plot(pg_plot)

  iter_plot = ggplot() +
    geom_point(aes( y= pg_out$relative_iterate_history, x = 1:length(pg_out$duality_gap_history))) +
    ggtitle(paste("PG Relative Iterate History: ", expression(lambda), " =", lambda , ", k = ", k )) +
    ylab("Relative Iterate") +
    xlab("Iteration")

  obj_plot = ggplot() +
    geom_point(aes( y= pg_out$relative_objective_history, x = 1:length(pg_out$relative_objective_history))) +
    ggtitle(paste("PG Relative Objective History: ", expression(lambda), " =", lambda , ", k = ", k )) +
    ylab("Relative Objective") +
    xlab("Iteration")

  dg_plot = ggplot() +
    geom_point(aes( y= pg_out$duality_gap_history, x = 1:length(pg_out$duality_gap_history))) +
    ggtitle(paste("PG Duality Gap History: ", expression(lambda), " =", lambda , ", k = ", k )) +
    ylab("Duality Gap") +
    xlab("Iteration")

  kkt_plot = ggplot() +
    geom_point(aes( y= pg_out$kkt_residual_history, x = 1:length(pg_out$duality_gap_history))) +
    ggtitle(paste("PG KKT Residual History: ", expression(lambda), " =", lambda , ", k = ", k )) +
    ylab("KKT Residual") +
    xlab("Iteration")

  multiplot(iter_plot, obj_plot, dg_plot, kkt_plot, cols=2)
}

apple_data = get.hist.quote( instrument = "AAPL",
                             start = "2018-01-01",
                             end = "2020-01-01",
                             quote="Close")

times = time(apple_data)
time_stock = cbind(times, as.data.frame(apple_data))

y=time_stock[,2]


k = 2
lambda = 0.1
max_iter = 1e2
tol=1e-3

make_plots(y=y,k=k, lambda = lambda, max_iter=max_iter, tol=tol)


