\name{dcElerian}
\alias{dcElerian}
\title{Approximated conditional law of a diffusion process by Elerian's method}
\description{
Approximated conditional densities for \code{X(t) | X(t0) = x0} of a diffusion process}
\usage{
dcElerian(x, t, x0, t0, theta, d, s, sx, log=FALSE)
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
  \item{sx}{partial derivative wrt \code{x} of the
      diffusion coefficient. See details.}
}
\details{
This function returns the value of the conditional density of
\code{X(t) | X(t0) = x0} at point \code{x}. 

All the functions \code{d}, \code{s}
and \code{sx} must be functions of \code{t}, \code{x} and \code{theta}.
}
\value{
  \item{x}{a numeric vector}
}
\author{Stefano Maria Iacus}
\references{
Elerian, O. (1998) A note on the existence of a closed form conditional 
density for the Milstein scheme, \emph{Working Paper, Nuffield College, 
Oxford University}. Available at \url{http://www.nuff.ox.ac.uk/economics/papers/}
}
\keyword{datagen}
\keyword{ts}
