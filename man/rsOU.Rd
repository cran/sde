\name{rsOU}
\alias{rsOU}
\alias{dsOU}
\alias{psOU}
\alias{qsOU}
\title{Ornstein-Uhlenbeck or Vasicek process stationary law}
\description{Density, distribution function, quantile function, and 
random generation for the stationary law of the Ornstein-Uhlenbeck process
also known as the Vasicek process.}
\usage{
dsOU(x, theta, log = FALSE)
psOU(x, theta, lower.tail = TRUE, log.p = FALSE) 
qsOU(p, theta, lower.tail = TRUE, log.p = FALSE)
rsOU(n=1, theta)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{theta}{parameter of the Ornstein-Uhlenbeck process; see details.}
  \item{n}{number of random numbers to generate from the conditional distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities \eqn{p}{p} are given as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \code{P[X <= x]}; 
  otherwise \code{P[X > x]}.}
}
\details{
This function returns quantities related to the stationary law
of the process solution of
\deqn{{\rm d}X_t = (\theta_1-\theta_2 X_t){\rm d}t + \theta_3 {\rm d}W_t.}{dX_t = (theta[1]-theta[2]*Xt)*dt + theta[3]*dWt.}

Contraints: \eqn{theta_2>0, \theta_3>0}{theta[2]>0, theta[3]>0}.

Please note that the process is stationary only if \eqn{\theta_2>0}{theta[2]>0}.
}
\value{
  \item{x}{a numeric vector}
}
\references{Uhlenbeck, G. E.,  Ornstein, L. S. (1930) On the theory of Brownian motion, 
\emph{Phys. Rev.}, 36, 823-841.

Vasicek, O. (1977) An Equilibrium Characterization of the Term 
Structure, \emph{Journal of Financial Economics},  5, 177-188. }
\author{Stefano Maria Iacus}
\seealso{\code{\link{rcOU}}}
\examples{
rsOU(n=1, theta=c(0,2,1))
}
\keyword{datagen}
\keyword{ts}
