\name{dcEuler}
\alias{dcEuler}
\title{Approximated conditional law of a diffusion process}
\description{
Approximated conditional densities for X(t)|X(t0)=x0 of a diffusion process}
\usage{
dcEuler(x, t, x0, t0, theta, d, s, log=FALSE)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{t}{lag or time.}
  \item{x0}{the value of the process at time \code{t0}. See details.}
  \item{t0}{intial time.}
  \item{theta}{parameter of the process. See details.}
  \item{log}{logical; if TRUE, probabilities p are given as log(p).}
  \item{d}{drift coefficient as a function. See details.}
  \item{s}{diffusion coefficient as a function. See details.}
}
\details{
This function returns the value of the conditional density of
\code{X(t) | X(t0) = x0} at point \code{x}. 

All the functions \code{d}, \code{s}, must be functions of \code{t0}, 
\code{x0} and \code{theta}.
}
\value{
  \item{x}{a numeric vector}
}
\author{Stefano Maria Iacus}
\note{This package is a companion to the book \emph{Simulation and Inference
for Stochastic Differential Equation}, Springer, NY.
}
\keyword{datagen}
\keyword{ts}
