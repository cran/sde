\name{rcBS}
\alias{rcBS}
\alias{dcBS}
\alias{pcBS}
\alias{qcBS}
\title{Black-Scholes-Merton or geometric Brownian motion process conditional law}
\description{Density, distribution function, quantile function, and 
random generation for the conditional law \eqn{X(t) | X(0) = x_0}{X(t) | X(0) = x0} 
of the Black-Scholes-Merton process
also known as the geometric Brownian motion process.}
\usage{
dcBS(x, Dt, x0, theta, log = FALSE)
pcBS(x, Dt, x0, theta, lower.tail = TRUE, log.p = FALSE) 
qcBS(p, Dt, x0, theta, lower.tail = TRUE, log.p = FALSE)
rcBS(n=1, Dt, x0, theta)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{Dt}{lag or time.}
  \item{x0}{the value of the process at time \code{t}; see details.}
  \item{theta}{parameter of the Black-Scholes-Merton process; see details.}
  \item{n}{number of random numbers to generate from the conditional distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities \eqn{p}{p} are given as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \code{P[X <= x]}; 
  otherwise, \code{P[X > x]}.}
}
\details{
This function returns quantities related to the conditional law
of the process solution of
\deqn{{\rm d}X_t = \theta_1 X_t {\rm d}t + \theta_2 X_t {\rm d}W_t.}{dX_t = theta[1]*Xt*dt + theta[2]*Xt*dWt.}

Constraints: \eqn{\theta_3>0}{theta[3]>0}.

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
\examples{
rcBS(n=1, Dt=0.1, x0=1, theta=c(2,1))
}
\keyword{datagen}
\keyword{ts}
