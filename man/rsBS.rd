\name{rsBS}
\alias{rsBS}
\alias{dsBS}
\alias{psBS}
\alias{qsBS}
\title{Black-Scholes-Merton or Geometric Brownian Motion process stationary law}
\description{Density, distribution function, quantile function and 
random generation for the stationary law of for the Blac-Scholes process
also known as Geometric Brownian Motion process}
\usage{
dsBS(x, theta, log = FALSE)
psBS(x, theta, lower.tail = TRUE, log.p = FALSE) 
qsBS(p, theta, lower.tail = TRUE, log.p = FALSE)
rsBS(n=1, theta)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{theta}{parameter of the Black-Scholes-Merton process. See details.}
  \item{n}{number of random numbers to generate from the conditional distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x], 
  otherwise, P[X > x].}
}
\details{
This function returns quantities related to the stationary law
of the process solution of
\code{dX_t = theta[1]*Xt*dt + theta[2]*Xt*dWt}.

Constraints: \code{theta[3]>0}.

}
\value{
  \item{x}{a numeric vector}
}
\references{ Black, F.,  Scholes, M.S. (1973) The pricing of options 
and corporate liabilities, \emph{Journal of Political Economy}, 81, 637-654.

Merton, R. C. (1973) Theory of rational option pricing, 
\emph{Bell Journal of Economics and Management Science}, 4(1), 141-183.
}
\author{Stefano Maria Iacus}
\note{This package is a companion to the book `Simulation and Inference
for Stochastic Differential Equation, Springer, NY.
}
\seealso{\code{\link{rcBS}}}
\examples{
rsBS(n=1, theta=c(2,1))
}
\keyword{datagen}
\keyword{ts}
