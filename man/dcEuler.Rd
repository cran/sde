\name{dcEuler}
\alias{dcEuler}
\title{Approximated conditional law of a diffusion process}
\description{
Approximated conditional densities for \eqn{X(t) | X(t_0) = x_0}{X(t) | X(t0) = x0} of a diffusion process.}
\usage{
dcEuler(x, t, x0, t0, theta, d, s, log=FALSE)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{t}{lag or time.}
  \item{x0}{the value of the process at time \code{t0}; see details.}
  \item{t0}{initial time.}
  \item{theta}{parameter of the process; see details.}
  \item{log}{logical; if TRUE, probabilities \eqn{p}{p} are given as \eqn{\log(p)}{log(p)}.}
  \item{d}{drift coefficient as a function; see details.}
  \item{s}{diffusion coefficient as a function; see details.}
}
\details{
This function returns the value of the conditional density of
\eqn{X(t) | X(t_0) = x_0}{X(t) | X(t0) = x0} at point \code{x}. 

The functions \code{d} and \code{s} must be functions of \code{t}, 
\code{x}, and \code{theta}.
}
\value{
  \item{x}{a numeric vector}
}
\author{Stefano Maria Iacus}
\keyword{datagen}
\keyword{ts}
