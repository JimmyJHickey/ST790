#' Compute the kth order trend filtering estimator
#'
#' @param y noisy signal
#' @param k order of differencing operator
#' @param lambda regularization parameter
#' @export
myTrendFilter <- function(y, k, lambda) {
 library(gurobi)

 n = length(y)

 Dk = myGetDkn(k = k, n = n)

 nrow_d = nrow(Dk)

 neg_id_d = -1 * Diagonal(nrow_d)
 pos_id_d = Diagonal(nrow_d)
 half_id_n = 0.5 * Diagonal(n)

 zero_d_n = Matrix(0, nrow=nrow_d, ncol = n, sparse=TRUE)
 zero_d_d = Matrix(0, nrow=nrow_d, ncol = nrow_d, sparse=TRUE)
 zero_n_d = Matrix(0, nrow=n, ncol = nrow_d, sparse=TRUE)
 zero_n_n = Matrix(0, nrow=n, ncol = n, sparse=TRUE)

 lambda_d = rep(lambda, nrow_d)

 linear = c(-1 * t(y), lambda_d, lambda_d)

 quadratic = rbind(
                   cbind(half_id_n, zero_n_d, zero_n_d),
                   cbind(zero_d_n, zero_d_d, zero_d_d),
                   cbind(zero_d_n, zero_d_d, zero_d_d)
                  )

 equality_constraints = cbind(Dk, neg_id_d, pos_id_d)
 inequality_constraints = rbind(
                               cbind(zero_d_n, neg_id_d, zero_d_d),
                               cbind(zero_d_n, zero_d_d, neg_id_d)
                               )

 model = list()


 model$A = rbind(equality_constraints, inequality_constraints)
 model$Q = quadratic
 model$obj = linear
 model$rhs = c(
               rep(0, nrow_d),
               rep(0, nrow_d),
               rep(0, nrow_d)
              )
 model$sense = c(
                 rep("=", nrow_d),
                 rep("<=", nrow_d),
                 rep("<=", nrow_d)
                )

 cat(dim( model$A))
 cat(" ")
 cat(dim( model$Q))

 result = gurobi(model)

 return(result)

}
