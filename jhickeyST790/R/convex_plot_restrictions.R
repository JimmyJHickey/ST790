###
#
# Plot restrictions to a line of a convex function.
#
# Jimmy Hickey
# 2020-01-22
#
###



#' @param fx handle to function in question
#' @param rx handle to function that returns a sequence of feasible sequence t and x + tv
#' @param nRow number of row plots
#' @param nCol number of column plots
plot_restrictions <- function(fx, rx, nRow=3, nCol=3) {
  par(mfrow=c(nRow, nCol))
  for (i in 1:(nRow*nCol)) {
    rand_rest <- rx()
    t <- rand_rest$t
    Z <- rand_rest$Z
    nt <- length(t)
    g <- double(nt)
    for (i in 1:nt) {
      g[i] <- fx(Z[[i]])
    }
    plot(t,g,type='l')
  }
}



#' Create restrictions of positive definite matrices.
#'
#' @param n Number of columns and rows
#' @param nt Number of restrictions to make
#'
#' @return
restrictions_positive_definite_matrices <- function(n = 5, nt = 1e2) {
  ## Create positive definite X
  svdS <- svd(matrix(rnorm(n**2),n,n))
  U <- svdS$u
  X <- U %*% diag(1+rnorm(n)**2) %*% t(U)

  ## Create positive definite V
  svdS <- svd(matrix(rnorm(n**2),n,n))
  U <- svdS$u
  V <- U %*% diag(1+rnorm(n)**2) %*% t(U)

  ## Create sequence positive increasing sequence t and Z
  Z <- vector(mode="list", length=nt)
  t <- cumsum(runif(nt))
  for (i in 1:nt) {
    Z[[i]] <- X + t[i]*V
  }
  return(list(t=t, Z=Z))
}
