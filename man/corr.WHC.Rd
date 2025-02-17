% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corr.WHC.R
\name{corr.WHC}
\alias{corr.WHC}
\title{Correlation adjusted weighted Hochberg method}
\usage{
corr.WHC(p, w, corr.mat, a = 0.5, b = 0.6, penalty = NULL, alpha = 0.05)
}
\arguments{
\item{p}{A numeric vector. A length \eqn{m} P-value vector from multiple tests.}

\item{w}{A numeric vector. Any non-negative real numbers to denote the
importance of the endpoints. Length must be equal to \eqn{m}. A single value,
e.g. \code{w = 1}, represents equal weight. \code{WHC} can scale the weight vector as if
the sum of weight is not 1.}

\item{corr.mat}{A matrix. The dimension must be \eqn{m \times m}. Positive correlation
is the theoretical assumption, however, it is robust to run with some negative
elements in the correlation matrix.}

\item{a}{A numeric number. \eqn{a \in (0,1)} determines the greatest penalty on weight,
Default \code{a=0.5}. Details see Huang et al (2024+).}

\item{b}{A numeric number. \eqn{b \in (0,1)} is the shape parameter of the penalty function.
\eqn{b = 1} produce a linear function.}

\item{penalty}{A function. User can define their own penalty function.
The basic rule is the function must be monotone decreasing from 0 to 1,
and range from 1 to \eqn{a} where \eqn{a \in (0,1)}. A convex function is recommended.
Concave function can produce result, but have no meaning on alpha conserving.}

\item{alpha}{A real number. \eqn{1-\alpha} is the confidence level, \code{alpha} must between (0, 1).}
}
\value{
A table contains p-values, weights, adjusted critical values, significance
}
\description{
A new method implement correlation correction based on weighted Hochberg. An ACF is
applied for weight reduction to conserve alpha. Details see Huang et al. (2024+).
A correlation structure with too many zero leads the method reduce to weighted
Hochberg.
}
\examples{
m <- 5
corr.WHC(
  p = runif(m),
  w = runif(m),
  corr.mat = cor(matrix(runif(10*m), ncol = m))
)
}
\references{
Huang, X. -W., Hua, J., Banerjee, B., Wang, X., Li, Q. (2024+). Correlated weighted Hochberg procedure. In-preparation.
}
