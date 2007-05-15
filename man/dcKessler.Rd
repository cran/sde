\name{dcKessler}
\alias{dcKessler}
\title{Approximated conditional law of a diffusion process by Kessler's method}
\description{
Approximated conditional densities for X(t)|X(t0)=x0 of a diffusion process}
\usage{
dcKessler(x, t, x0, t0, theta, d, dx, dxx, s, sx, sxx, log=FALSE) 
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
  \item{s}{diffusion coefficient as a function. See details.}
  \item{sx}{partial derivative wrt \code{x} of the
      diffusion coefficient. See details.}
  \item{sxx}{second partial derivative wrt \code{x^2} of the
      diffusion coefficient. See details.}
}
\details{
This function returns the value of the conditional density of
\code{X(t) | X(t0) = x0} at point \code{x}. 

All the functions \code{d}, \code{dx}, \code{dxx}, \code{dt}, \code{s}, \code{sx} 
and \code{sxx} must be functions of \code{t0}, \code{x0} and \code{theta}.
}
\value{
  \item{x}{a numeric vector}
}
\author{Stefano Maria Iacus}
\references{
Kessler, M. (1997)  Estimation of an ergodic diffusion from discrete observations, 
\emph{Scand. J. Statist.}, 24, 211-229.
}
\note{This package is a companion to the book \emph{Simulation and Inference
for Stochastic Differential Equation}, Springer, NY.
}
\keyword{datagen}
\keyword{ts}
