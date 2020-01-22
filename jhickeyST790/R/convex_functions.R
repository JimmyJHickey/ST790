#' Take the log of the determinant of a matrix.
#'
#' @param in_matrix A square, positive definite matrix.
#' @return The log of the determinant of \code{in_matrix}.
#' @examples
logdet <- function(in_matrix) {

  # check that the input matrix is square
  if (nrow(in_matrix) != ncol(in_matrix))
    stop("'X' must be a square matrix")

  matrix_det <- det(in_matrix)

  if (matrix_det <= 0)
    stop("'in_matrix' must be positive definite")

  return(-log(matrix_det))
}

