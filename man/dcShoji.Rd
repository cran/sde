\name{dcShoji}
\alias{dcShoji}
\title{Approximated conditional law of a diffusion process by Shoji-Ozaki method}
\description{
Approximated conditional densities for \code{X(t) | X(t0) = x0} of a diffusion process}
\usage{
dcShoji(x, t, x0, t0, theta, d, dx, dxx, dt, s, log=FALSE)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{t}{lag or time.}
  \item{x0}{the value of the process at time \code{t0}. See details.}
  \item{t0}{intial time.}
  \item{theta}{parameter of the process. See details.}
  \item{log}{logical; if TRUE, probabilities p are given as log(p).}
  \item{d}{drift coefficient as a function. See details.}
  \item{dx}{partial derivative wrt \code{x} of the
      drift coefficient. See details.}
  \item{dxx}{second partial derivative wrt \code{x^2} of the
      drift coefficient. See details.}
  \item{dt}{partial derivative wrt \code{t} of the
      drift coefficient. See details.}
  \item{s}{diffusion coefficient as a function. See details.}
}
\details{
This function returns the value of the conditional density of
\code{X(t) | X(t0) = x0} at point \code{x}. 

All the functions \code{d}, \code{dx}, \code{dxx}, \code{dt} and \code{s} 
must be functions of \code{t}, \code{x} and \code{theta}.
}
\value{
  \item{x}{a numeric vector}
}
\author{Stefano Maria Iacus}
\references{
Shoji, L., Ozaki, T. (1998) Estimation for nonlinear stochastic differential
equations by a local linearization method, 
\emph{Stochastic Analysis and Applications}, 16, 733-752.
}
\keyword{datagen}
\keyword{ts}
