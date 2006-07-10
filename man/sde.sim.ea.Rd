\name{sde.sim.ea}
\alias{sde.sim.ea}
\title{Exact simulation algorithm EA1}
\description{
A rejection-method scheme for diffusions.
}
\usage{
sde.sim.ea(t0 = 0, T = 1, X0 = 1, N = 100, delta, drift, drift.x, k1, k2, phi, max.psi = 1000, rh, A)
}
\arguments{
  \item{t0}{time origin}
  \item{T}{horizon of simulation}
  \item{X0}{initial value of the process}
  \item{N}{number of simulation steps}
  \item{delta}{time-step of the simulation}
  \item{drift}{drift coeffcient: a expression of two variables \code{t} and \code{x}}
  \item{drift.x}{partial derivative of drift coeffcient wrt to \code{x}: a function of two variables \code{t} and \code{x}}
  \item{k1}{lower bound for \code{psi(x) = 0.5*drift(x)^2 + 0.5*drift.x(x)}}
  \item{k2}{upper bound for \code{psi(x) = 0.5*drift(x)^2 + 0.5*drift.x(x)}}
  \item{phi}{the function \code{psi(x) - k1}}
  \item{max.psi}{upper value of the support of \code{psi} to search for its maximum}
  \item{rh}{the rejection function. Seed details.}
  \item{A}{\code{A(x)} is the integral of the \code{drift} between \code{0} and \code{x}}
}
\details{
	If missing \code{k1} and \code{k2} as well as \code{A}, \code{rh} and \code{phi}
	are calculated numerically by the function.
}
\value{
  \item{X}{an invisible \code{ts} object}
}
\references{Beskos, A. and Roberts, G.O. (2005) Exact simulation of diffusions, \emph{Ann. Appl. Probab.},  4, 2422-2444.}
\author{Stefano Maria Iacus}
\seealso{\code{\link{sde.sim}}}
\examples{
set.seed(123)
d <- expression(sin(x))
d.x <- expression(cos(x)) 
A <- function(x) 1-cos(x)
sde.sim.ea(delta=1/20, X0=0, N=500, drift=d, drift.x = d.x, A=A) -> Y
plot(Y, main="Periodic drift")
}
\keyword{datagen}
\keyword{ts}
