#' Weighted Hochberg method
#'
#' `WHC` performs the weighted Hochberg method proposed by Tamhane and Liu (2008).
#'
#' @param p A numeric vector. A length \eqn{m} P-value vector from multiple tests.
#' @param w A numeric vector. Any non-negative real numbers to denote the
#' importance of the endpoints. Length must be equal to \eqn{m}. A single value,
#' e.g. `w = 1`, represents equal weight. `WHC` can scale the weight vector as if
#' the sum of weight is not 1.
#' @param alpha A real number. \eqn{1-\alpha} is the confidence level, `alpha` must between (0, 1).
#' @returns A table contains p-values, weights, adjusted critical values, significance
#' @references Tamhane, A. C., & Liu, L. (2008). On weighted Hochberg procedures. Biometrika, 95(2), 279-294.
#' @examples
#' m <- 5
#' WHC(
#'   p = runif(m),
#'   w = runif(m)
#' )
#'@export
WHC <- function(
    p,
    w,
    alpha = 0.05
){

  m <- length(p)

  if (length(w)==1) {
    w0 <- rep(1/m, m)
  } else {
    w0 <- w/sum(w)
  }

  out <- tibble::tibble(p = p, w = w0) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    dplyr::arrange(dplyr::desc(p)) %>%
    dplyr::mutate(
      w.alpha = w / (cumsum(w)) * alpha,
      sig = as.numeric(p<=w.alpha)
    ) %>%
    dplyr::arrange(id) %>%
    dplyr::rename(
      "adj.alpha" = "w.alpha"
    )

  return(
    out %>% dplyr::select(p, w, adj.alpha, sig)
  )

}
