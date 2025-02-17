% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corrmat_block.R
\name{corrmat_block}
\alias{corrmat_block}
\title{Block design correlation matrix}
\usage{
corrmat_block(d, rho)
}
\arguments{
\item{d}{An integer vector. Length \eqn{B} of block dimensions. Element of \code{d}
can be \code{1}, it would not generate a sub-matrix with the corresponding element in
\code{rho}, but just \eqn{1}.}

\item{rho}{A numeric vector. A length \eqn{B} vector of correlation coefficients,
represent \eqn{B} different block of correlation matrix.}
}
\value{
A correlation matrix
}
\description{
An easy function to generate a block design correlation matrix. Each diagonal
element \eqn{R_i} is a compound symmetric matrix with dimension
\eqn{d_i \times d_i}. Correlation coefficient in each block is \eqn{\rho_i}.
All the off-diagonal elements are \eqn{0}.
}
\examples{
corrmat_block(
  d = c(2,3,4),
  rho = c(0.1, 0.3, 0.5)
)
}
