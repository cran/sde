\name{simple.ef}
\alias{simple.ef}
\title{Simple estimating functions of types I and II}
\description{Apply a simple estimating function to find estimates of the
parameters of a process solution of a stochastic differential equation.
}
\usage{
simple.ef(X, f, guess, lower, upper) 
}
\arguments{
  \item{X}{a ts object containing a sample path of an sde.}
  \item{f}{a list of expressions of \code{x} and/or \code{y} and the parameters to be estimated; see details.}
  \item{guess}{initial value of the parameters; see details.}
  \item{lower}{lower bounds for the parameters; see details.}
  \item{upper}{upper bounds for the parameters; see details.}
}
\details{
The function \code{simple.ef} minimizes a simple estimating function
of the form \code{sum_i f_i(x,y;theta) = 0} or \code{sum_i f_i(x;theta)}
as a function of \code{theta}. The index \code{i} varies in \code{1:length(theta)}.

The list \code{f} is a list of expressions in \code{x} or \code{(x,y)}.
}
\value{
  \item{x}{a vector of estimates}
}
\references{
Kessler, M. (1997)  Estimation of an ergodic diffusion from discrete observations,
 \emph{Scand. J. Statist.}, 24, 211-229.

Kessler, M. (2000) Simple and Explicit Estimating Functions for a Discretely Observed 
Diffusion Process, \emph{Scand. J. Statist.}, 27, 65-82. 
}
\author{Stefano Maria Iacus}
\examples{
set.seed(123); 
# Kessler's estimator for O-H process
K.est <- function(x) {
  n.obs <- length(x)
  n.obs/(2*(sum(x^2)))
}

# Least squares estimators for the O-H process
LS.est <- function(x) {
  n <- length(x) -1
  k.sum <- sum(x[1:n]*x[2:(n+1)])
  dt <- deltat(x)
  ifelse(k.sum>0, -log(k.sum/sum(x[1:n]^2))/dt, NA)
}

d <- expression(-1 * x)
s <- expression(1) 
x0 <- rnorm(1,sd=sqrt(1/2))
sde.sim(X0=x0,drift=d, sigma=s,N=2500,delta=0.1) -> X
 
# Kessler's estimator as estimating function
f <- list(expression(2*theta*x^2-1))
simple.ef(X, f, lower=0, upper=Inf)
K.est(X)

# Least Squares estimator as estimating function
f <- list(expression(x*(y-x*exp(-0.1*theta))))
simple.ef(X, f, lower=0, upper=Inf)
LS.est(X)
}
\keyword{ts}
