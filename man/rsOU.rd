\name{rsOU}
\alias{rsOU}
\alias{dsOU}
\alias{psOU}
\alias{qsOU}
\title{Ornstein-Uhlenbeck or Vasicek process stationary law}
\description{Density, distribution function, quantile function and 
random generation for the stationary law of for the Ornstein-Uhlenbeck process
also known as Vasicek process}
\usage{
dsOU(x, theta, log = FALSE)
psOU(x, theta, lower.tail = TRUE, log.p = FALSE) 
qsOU(p, theta, lower.tail = TRUE, log.p = FALSE)
rsOU(n=1, theta)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{theta}{parameter of the Ornstein-Uhlenbeck process. See details.}
  \item{n}{number of random numbers to generate from the conditional distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities p are given as log(p).}
  \item{lower.tail}{logical; if TRUE (default), probabilities are P[X <= x], 
  otherwise, P[X > x].}
}
\details{
This function returns quantities related to the stationary law
of the process solution of
\code{dX_t = theta[1]*(theta[2] - Xt)*dt + theta[3]*dWt}.

Contraints: \code{theta[1]>0, theta[3]>0}.

Please note that the process is stationary only if \code{theta[1]>0}.
}
\value{
  \item{x}{a numeric vector}
}
\references{Uhlenbeck, G. E.,  Ornstein, L. S. (1930) On the theory of Brownian motion, 
\emph{Phys. Rev.}, 36, 823-841.

Vasicek, O. (1977) An Equilibrium Characterization of the Term 
Structure, \emph{Journal of Financial Economics},  5, 177-188. }
\author{Stefano Maria Iacus}
\note{This package is a companion to the book `Simulation and Inference
for Stochastic Differential Equation, Springer, NY.
}
\seealso{\code{\link{rcOU}}}
\examples{
rsOU(n=1, theta=c(2,0,1))
}
\keyword{datagen}
\keyword{ts}
