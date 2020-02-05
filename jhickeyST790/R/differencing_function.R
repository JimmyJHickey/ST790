#' Compute kth order differencing matrix
#'
#' @param k order of the differencing matrix
#' @param n Number of time points
#' @return kth order differencing matrix
myGetDkn <- function(k, n)
{
  # create D1n
  if (k==1)
  {
    D1 = matrix(nrow = n-1, ncol = n)
    zeros_matrix = rep(0, n)

    for(i in 1:n-1)
    {
      D1[i,] = zeros_matrix
      D1[i,i] = -1
      D1[i, i+1] = 1
    } # for
    return(D1)
  } # endif

  return(myGetDkn(k=1, n=n - k + 1) %*% myGetDkn(k=k-1, n=n))

}

