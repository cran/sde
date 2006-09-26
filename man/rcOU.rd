\name{rcBS}
\alias{rcBS}
\alias{dcBS}
\alias{pcBS}
\alias{qcBS}
\title{Black-Scholes-Merton or Geometric Brownian Motion process conditional law}
\description{Density, distribution function, quantile function and 
random generation for the conditional law Xt|X0=x0 of the Black-Scholes-Merton process
also known as Geometric Brownian Motion process}
\usage{
dcBS(x, t, x0, theta, log = FALSE)
pcBS(x, t, x0, theta, lower.tail = TRUE, log.p = FALSE) 
qcBS(p, t, x0, theta, lower.tail = TRUE, log.p = FALSE)
rcBS(n=1, t, x0, theta)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{t}{lag or time}
  \item{x0}{the value of the process at time \code{t=0}. See details.}
  \item{theta}{parameter of the Black-Scholes-Merton process. See details.}
  \item{n}{number of random numbers to generate from the conditional distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x], 
  otherwise, P[X > x].}
}
\details{
This function returns quantities related to the conditional law
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
\note{This package is a companion to the book \emph{Simulation and Inference
for Stochastic Differential Equation}, Springer, NY.
}
\seealso{\code{\link{rsBS}}}
\examples{
rcBS(n=1, t=0.1, x0=1, theta=c(2,1))
}
\keyword{datagen}
\keyword{ts}
