\name{simple.ef2}
\alias{simple.ef2}
\title{Simple estimating function based on the infinitesimal generator a the diffusion process}
\description{Apply a simple estimating function based on the infinitesimal
generator of a diffusion to find estimates of the
parameters of a process solution of that particular stochastic differential equation.
}
\usage{
simple.ef2(X, drift, sigma, h, h.x, h.xx, guess, lower, 
           upper) 
}
\arguments{
  \item{X}{a \code{ts} object containing a sample path of an sde.}
  \item{drift}{an expression for the drift coefficient; see details.}
  \item{sigma}{an expression for the diffusion coefficient; see details.}
  \item{h}{an expression of \code{x} and the parameters to be estimated; see details.}
  \item{h.x}{an expression of \code{x} containing the first derivative of \code{h}; see details.}
  \item{h.xx}{an expression of \code{x} containing the second derivative of \code{h}; see details.}
  \item{guess}{initial value of the parameters; see details.}
  \item{lower}{lower bounds for the parameters; see details.}
  \item{upper}{upper bounds for the parameters; see details.}
}
\details{
The function \code{simple.ef2} minimizes the simple estimating function
of the form \code{sum_i f_i(x;theta) = 0}, where \code{f} is the result of
 applying  the infinitesimal generator of the diffusion to the
function \code{h}. This involves the drift and diffusion coefficients plus
the first two derivatives of \code{h}. If not provided by the user, the derivatives
are calculated by the function.

}
\references{
Kessler, M. (1997)  Estimation of an ergodic diffusion from discrete observations,
 \emph{Scand. J. Statist.}, 24, 211-229.

Kessler, M. (2000) Simple and Explicit Estimating Functions for a Discretely Observed 
Diffusion Process, \emph{Scand. J. Statist.}, 27, 65-82. 
}
\value{
  \item{x}{a vector of estimates}
}
\author{Stefano Maria Iacus}
\examples{
set.seed(123) 
d <- expression(10 - x)
s <- expression(sqrt(x)) 
x0 <- 10
sde.sim(X0=x0,drift=d, sigma=s,N=1500,delta=0.1) -> X

# rather difficult problem unless a good initial guess is given
d <- expression(alpha + theta*x)
s <- expression(x^gamma) 
h <- list(expression(x), expression(x^2), expression(x^2))
simple.ef2(X, d, s, h, lower=c(0,-Inf,0), upper=c(Inf,0,1))
}
\keyword{ts}
