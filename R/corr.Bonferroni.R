#' ICC adjusted Bonferroni method
#'
#' `corr.Bonferroni` performs the ICC adjusted Bonferroni method proposed by
#' Shi, Pavey, and Carter(2012). Power law approximation by `r` is tricky, suggested
#' options was listed in the paper.
#'
#' @param p A numeric vector. A length \eqn{m} P-value vector from multiple tests.
#' @param ICC A number. Intraclass correlation correction factor, a real number between (0, 1).
#' @param r A number. Tuning parameter for g** between (0, 1). Default `r=0`.
#' @param alpha A real number. \eqn{1-\alpha} is the confidence level, `alpha` must between (0, 1).
#' @returns A numeric vector of adjusted p-values.
#' @references Shi, Q., Pavey, E. S., & Carter, R. E. (2012). Bonferroni‐based correction factor for multiple, correlated endpoints. Pharmaceutical statistics, 11(4), 300-309.
#' @examples
#' m <- 10
#' corr.Bonferroni(
#'   p = runif(m),
#'   ICC = 0.3
#' )
#'@export

corr.Bonferroni <- function(
    p,
    ICC,
    r = 0,
    alpha = 0.05
) {
  g <- length(p)

  g_1 <- (g + 1) - (1 + (g - 1)*ICC) # basic adjustment

  g_2 <- g_1^( (1 + r)^(g-g_1) ) # advanced adjustment g_1 < g_2 < g

  return(
    pmin(g_2 * p, rep(1,g))
  )

}
