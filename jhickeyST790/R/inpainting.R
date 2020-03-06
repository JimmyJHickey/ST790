#' Create an n x n sparse diagonal matrix
#'
#' @param n size of matrix
#' @param diag number on diag
#' @export
sparseDiag = function(n, diag)
{
 return(bandSparse(n, n, 0, as.matrix(rep(diag,n))))
}


#' Create an m x n sparse 0  matrix
#'
#' @param nrow number of rows
#' @param ncol number of columns
#' @export?
sparseZero = function(nrow, ncol)
{
 return(bandSparse(nrow, ncol, 0, as.matrix(rep(0,nrow))))
}

#' Compute the inpainting differencing matrix
#'
#' @param m number of rows
#' @param n number of columns
#' @export
myGetTV2d <- function(m,n) {

 In = sparseDiag(n, 1)
 Im = sparseDiag(m, 1)

 return( rbind(kronecker(In, myGetDkn(1, m, negative=TRUE)),
               kronecker(myGetDkn(1, n, negative=TRUE),Im)))
}

#' Remove missing data to calculate y and P_Omega
#'
#' @param inMatrix input data matrix
#' @export
removeMissing = function(inMatrix)
{
 in_vec = c(inMatrix)
 mn = length(in_vec)

 not_na_indices = which(!is.na(inMatrix) & in_vec != -1 )
 num_not_na = length(not_na_indices)

 y = in_vec[not_na_indices]
 P_Omega = sparseMatrix(i = 1:num_not_na, j = not_na_indices, x = 1, dims=c(num_not_na, mn))

 out_list = list(
                 "p_omega" = P_Omega,
                 "y" = y
              )

 return(out_list)
}


#' Compute the inpainting completion
#'
#' @param Y data matrix
#' @export
myInpaint <- function(Y) {
 # x vector = (x, w+, w-)

 library(gurobi)

 m = nrow(Y)
 n = ncol(Y)
 mn = m * n


 Dhat = myGetTV2d(m, n)
 Dhat_nrow = nrow(Dhat)

 no_missing = removeMissing(Y)
 p_omega = no_missing$p_omega
 y = no_missing$y

 neg_Imn = sparseDiag(mn, -1)

 neg_Idrow = sparseDiag(Dhat_nrow, -1)
 zero_drow = sparseDiag(Dhat_nrow, 0)
 one_vec_drow =  rep(1, Dhat_nrow)


 zero_p_drow = sparseZero(nrow(p_omega), Dhat_nrow)
 zero_mn_drow = sparseZero(mn, Dhat_nrow)
 zero_drow_mn = sparseZero(Dhat_nrow, mn)
 zero_drow_drow = sparseZero(Dhat_nrow, Dhat_nrow)

 zero_vec_mn = rep(0, mn)
 zero_vec_drow = rep(0, Dhat_nrow)

 equality_constraints = rbind(
                           cbind(no_missing$p_omega, zero_p_drow, zero_p_drow),
                           cbind(Dhat, neg_Idrow, -1 * neg_Idrow)
                          )


 inequality_constraints = rbind(
                           cbind(neg_Imn, zero_mn_drow, zero_mn_drow),
                           cbind(zero_drow_mn, neg_Idrow, zero_drow_drow),
                           cbind(zero_drow_mn, zero_drow_drow, neg_Idrow)
                           )

 gurobi_model = list()

 gurobi_model$A = rbind(equality_constraints, inequality_constraints)
 gurobi_model$obj = c(zero_vec_mn, one_vec_drow, one_vec_drow)
 gurobi_model$modelsense = 'min'
 gurobi_model$rhs = c(y, zero_vec_drow,
                      zero_vec_mn, zero_vec_drow, zero_vec_drow)
 gurobi_model$sense = c(
                        rep('=', length(y)),
                        rep('=', length(zero_vec_drow)),
                        rep('<=', length(zero_vec_mn)),
                        rep('<=', length(zero_vec_drow)),
                        rep('<=', length(zero_vec_drow)))

 params <- list(OutputFlag=0)

 cat(length(gurobi_model$rhs))

 result <- gurobi(gurobi_model, params)

 return(result)
}
