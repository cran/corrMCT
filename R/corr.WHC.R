#' Correlation adjusted weighted Hochberg method
#'
#' A new method implement correlation correction based on weighted Hochberg. An ACF is
#' applied for weight reduction to conserve alpha. Details see Huang et al. (2024+).
#' A correlation structure with too many zero leads the method reduce to weighted
#' Hochberg.
#'
#'
#' @param p A numeric vector. A length \eqn{m} P-value vector from multiple tests.
#' @param w A numeric vector. Any non-negative real numbers to denote the
#' importance of the endpoints. Length must be equal to \eqn{m}. A single value,
#' e.g. `w = 1`, represents equal weight. `WHC` can scale the weight vector as if
#' the sum of weight is not 1.
#' @param corr.mat A matrix. The dimension must be \eqn{m \times m}. Positive correlation
#' is the theoretical assumption, however, it is robust to run with some negative
#' elements in the correlation matrix.
#' @param theta A numeric number. \eqn{\theta \in (0,1)} is the shape parameter of alpha conserving function.
#' \eqn{\theta = 0} is a straight line \eqn{f(\rho) = 1-\gamma\rho}. For \eqn{\theta > 0}, greater the value, more
#' aggressive decreasing it allows. Default `theta=0.10` is a relative conservative suggestion. Details see
#' Huang et al (2024+).
#' @param gamma A numeric number. \eqn{\gamma \in (0,1)} is the scale parameter of alpha conserving function.
#' \eqn{\gamma = 1} makes the function approach 1 when \eqn{\rho = 1}. \eqn{\gamma < 1} sets a bound to prevent
#' too much penalty on weight. Default `theta=0.50` is a moderate suggestion. Details see Huang et al (2024+).
#' @param ACF A function. User can define their own alpha conserving function.
#' The basic rule is the function must be monotone decreasing from 0 to 1,
#' and range from 1 to \eqn{a} where \eqn{a \in (0,1)}. A convex function is recommended.
#' Concave function can produce result, but have no meaning on alpha conserving.
#' @param alpha A real number. \eqn{1-\alpha} is the confidence level, `alpha` must between (0, 1).
#' @returns A table contains p-values, weights, adjusted critical values, significance
#' @references Huang, X. -W., Hua, J., Banerjee, B., Wang, X., Li, Q. (2024+). Correlated weighted Hochberg procedure. In-preparation.
#' @examples
#' m <- 5
#' corr.WHC(
#'   p = runif(m),
#'   w = runif(m),
#'   corr.mat = cor(matrix(runif(10*m), ncol = m))
#' )
#'@export

corr.WHC <- function(
    p,
    w,
    corr.mat,
    theta = 0.2,
    gamma = 0.3,
    ACF = NULL,
    alpha = 0.05
) {
  # browser()
  m <- length(p)

  if (length(w)==1) {
    w0 <- rep(1/m, m)
  } else {
    w0 <- w/sum(w)
  }

  temp <- tibble::tibble(p = p, w = w0) %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    dplyr::arrange(dplyr::desc(p)) %>%
    dplyr::mutate(adj.w = 0)

  ACF0 <- function(rho){
    (1-gamma*rho)^exp(theta*(m-2))
  }

  # ACF0 <- function(rho){
  #   (1-rho)^exp(1/4*log(2)*(m-2))
  # }

  if (is.null(ACF)) {
    ACF <- ACF0
  }

  for (i in 1:length(p)) {

    rho.order.vec <- c()
    for (k in 1:i) {
      # rho.order.vec[k] <- tryCatch(
      #   ifelse(i == k, 0, corr.mat[temp$id[i],temp$id[k]]),
      #   error = function(e) {browser()}
      # )
      rho.order.vec[k] <- ifelse(i == k, 0, corr.mat[temp$id[i],temp$id[k]])

    }

    temp$adj.w[i] <- temp$w[i] / sum(temp$w[1:i] * ACF(rev(rho.order.vec)))

  }

  out <- temp %>%
    dplyr::mutate(
      adj.alpha =  adj.w * alpha,
      sig = as.numeric(p<=adj.alpha)
    )

  for (i in 1:length(p)) {

    if (out$sig[i] == 1) {
      out$sig[i:length(p)] = 1
      break;
    }
  }

  out <- out %>% dplyr::arrange(id)

  return(
    out %>% dplyr::select(p, w, adj.alpha, sig)
  )

}
