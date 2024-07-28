#' Compound symmetric correlation matrix
#'
#' An easy function to generate a compound symmetric correlation matrix
#'
#' @param m An integer. Dimension of the correlation matrix.
#' @param rho A number. Correlation coefficient between \eqn{(0,1)}
#' @returns A correlation matrix
#' @examples
#' corrmat_CS(
#'   m = 3,
#'   rho = 0.2
#' )
#'@export
#'

corrmat_CS <- function(m, rho) {

  out <- corrmat_block(d = m[1], rho = rho[1])

  return(out)

}

