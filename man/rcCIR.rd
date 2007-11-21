\name{rcCIR}
\alias{dcCIR}
\alias{pcCIR}
\alias{qcCIR}
\alias{rcCIR}
\title{Conditional law of the Cox-Ingersoll-Ross process}
\description{
Density, distribution function, quantile function and 
random generation for the conditional law \eqn{X(t+D_t) | X(t)=x_0}{X(t+D_t) | X(t)=x0} of the Cox-Ingersoll-Ross
 process.}
\usage{
dcCIR(x, Dt, x0, theta, log = FALSE)
pcCIR(x, Dt, x0, theta, lower.tail = TRUE, log.p = FALSE) 
qcCIR(p, Dt, x0, theta, lower.tail = TRUE, log.p = FALSE)
rcCIR(n=1, Dt, x0, theta)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{Dt}{lag or time.}
  \item{x0}{the value of the process at time \code{t}; see details.}
  \item{theta}{parameter of the Ornstein-Uhlenbeck process; see details.}
  \item{n}{number of random numbers to generate from the conditional distribution.}
  \item{log, log.p}{logical; if TRUE, probabilities \eqn{p}{p} are given as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if TRUE (default), probabilities are \code{P[X <= x]}; 
  otherwise \code{P[X > x]}.}
}
\details{
This function returns quantities related to the conditional law
of the process solution of
\deqn{{\rm d}X_t = (\theta_1-\theta_2 X_t){\rm d}t + \theta_3\sqrt{X_t}{\rm d}W_t.}{dX_t = (theta[1]-theta[2]*Xt)*dt + theta[3]*sqrt(X_t)*dWt.}

Constraints: \eqn{2\theta_1> \theta_3^2}{2*theta[1]> theta[3]^2}, all \eqn{\theta}{theta} positive.
}
\value{
  \item{x}{a numeric vector}
}
\references{Cox, J.C., Ingersoll, J.E., Ross, S.A. (1985) A theory 
of the term structure of interest rates,  \emph{Econometrica}, 53, 385-408.}
\author{Stefano Maria Iacus}
\seealso{\code{\link{rsCIR}}}
\examples{
rcCIR(n=1, Dt=0.1, x0=1, theta=c(6,2,2))
}
\keyword{datagen}
\keyword{ts}
