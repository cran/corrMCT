#' Block AR(1) design correlation matrix
#'
#' An easy function to generate a block AR(1) design correlation matrix. Each diagonal
#' element \eqn{R_i} is an AR(1) correlation matrix with dimension
#' \eqn{d_i \times d_i}. Correlation coefficient in each block is \eqn{\rho_i}.
#' All the off-diagonal elements are \eqn{0}.
#'
#' @param d An integer vector. Length \eqn{B} of block dimensions. Element of `d`
#' can be `1`, it would not generate a sub-matrix with the corresponding element in
#' `rho`, but just \eqn{1}.
#' @param rho A numeric vector. A length \eqn{B} vector of correlation coefficients,
#' represent \eqn{B} different block of correlation matrix.
#'
#' @returns A correlation matrix
#' @examples
#' corrmat_blockAR1(
#'   d = c(2,3,4),
#'   rho = c(0.1, 0.3, 0.5)
#' )
#'@export
corrmat_blockAR1 <- function(d, rho) {

  mat <- list()
  for (i in seq_along(rho)){
    mat[[i]] <- corrMCT::corrmat_AR1(d[i], rho[i])
  }

  out <- Matrix::bdiag(mat) %>% as.matrix()

  return(out)

}
