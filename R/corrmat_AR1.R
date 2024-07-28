#' AR(1) correlation matrix
#'
#' An easy function to generate a AR(1) correlation matrix.
#'
#' @param m An integer. Dimension of the correlation matrix.
#' @param rho A number. Correlation coefficient between \eqn{(0,1)}
#' @returns A correlation matrix
#' @examples
#' corrmat_AR1(
#'   m = 3,
#'   rho = 0.2
#' )
#'@export


corrmat_AR1 <- function(m, rho) {

  temp <- abs(
    matrix(1:m[1] - 1, nrow = m[1], ncol = m[1], byrow = TRUE) - (1:m[1] - 1)
  )

  return(rho^temp)

}

