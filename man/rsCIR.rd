\name{rsCIR}
\alias{rsCIR}
\alias{dsCIR}
\alias{psCIR}
\alias{qsCIR}
\title{Cox-Ingersoll-Ross process stationary law}
\description{Density, distribution function, quantile function and 
random generation for the stationary law of for the Cox-Ingersoll-Ross process}
\usage{
dsCIR(x, theta, log = FALSE)
psCIR(x, theta, lower.tail = TRUE, log.p = FALSE) 
qsCIR(p, theta, lower.tail = TRUE, log.p = FALSE)
rsCIR(n=1, theta)
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
\code{dX_t = (theta[1] - theta[2]*Xt)*dt + theta[3]*sqrt(X_t)*dWt}.

Constraints: \code{2*theta[1] > theta[3]^2, theta's>0}.
}
\value{
  \item{x}{a numeric vector}
}
\references{Cox, J.C., Ingersoll, J.E., Ross, S.A. (1985) A theory 
of the term structure of interest rates,  \emph{Econometrica}, 53, 385-408.}
\author{Stefano Maria Iacus}
\note{This package is a companion to the book \emph{Simulation and Inference
for Stochastic Differential Equation}, Springer, NY.
}
\seealso{\code{\link{rsCIR}}}
\examples{
rsCIR(n=1, theta=c(6,2,1))
}
\keyword{datagen}
\keyword{ts}
