#' Compute kth order differencing matrix
#'
#' @param k order of the differencing matrix
#' @param n Number of time points
#' @param negative Flip signs, row one will look lie 1 -1 0 ... 0 (default: FALSE)
#' @return kth order differencing matrix
myGetDkn <- function(k, n, negative=FALSE)
{
  library(Matrix)
  neg = negative

  ii = 1 + (!negative) * (-2)
  iiplusone = -ii

   # create D1n
  if (k==1)
  {
    D1 = matrix(nrow = n-1, ncol = n)
    zeros_matrix = rep(0, n)

    for(i in 1:n-1)
    {
      D1[i,] = zeros_matrix
      D1[i,i] = ii
      D1[i, i+1] = iiplusone
    } # for
    return(Matrix(D1, sparse = TRUE))
  } # endif

  return(myGetDkn(k=1, n=n - k + 1, neg) %*% myGetDkn(k=k-1, n=n, neg))

}

