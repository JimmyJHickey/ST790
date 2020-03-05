#' Create an n x n sparse identity matrix
#'
#' @param n size of matrix
#' @export
sparseIdentity = function(n)
{
 return(bandSparse(n, n, 0, as.matrix(rep(1,n))))
}

#' Create an n x n sparse zero matrix
#'
#' @param n size of matrix
#' @export
sparseZero = function(n)
{
 return(bandSparse(n, n, 0, as.matrix(rep(0,n))))
}

#' Compute the inpainting differencing matrix
#'
#' @param m number of rows
#' @param n number of columns
#' @export
myGetTV2d <- function(m,n) {

 In = sparseIdentity(n)
 Im = sparseIdentity(m)

 return( rbind(kronecker(In, myGetDkn(1, m, negative=TRUE)),
               kronecker(myGetDkn(1, n, negative=TRUE),Im)))
}


#'
#' #' Compute the inpainting completion
#' #'
#' #' @param Y data matrix
#' #' @export
#' myInpaint <- function(Y) {
#'  # x vector = (x, w+, w-)
#'
#'  library(gurobi)
#'
#'  m = nrow(Y)
#'  n = ncol(Y)
#'
#'  Dhat = myGetTV2d(m, n)
#'  Dhat_nrow = nrow(Dhat)
#'
#'
#'
#'  Imn = sparseIdentity(mn)
#'  Idrow = sparseIdentity(Dhat_nrow)
#'
#'
#'  equality_constraints =
#'
#'
#'  gurobi_model = list()
#'
#'  gurobi_model$A = rbind(equality_constraints, inequality_constraints)
#'  gurobi_model$obj = c(0, 1, 1)
#'  gurobi_model$modelsense = 'min'
#'  gurobi_model$rhs = c()
#'  gurobi_model$sense = c('=', '=', '<', '<', '<')
#'
#'  params <- list(OutputFlag=0)
#'
#'  result <- gurobi(gurobi_model, params)
#'
#'  print('Solution:')
#'  print(result$objval)
#'  print(result$x)
#'
#'  # Clear space
#'  rm(model, result, params)
#'
#' }
