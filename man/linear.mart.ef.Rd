\name{linear.mart.ef}
\alias{linear.mart.ef}
\title{Linear martingale estimating function}
\description{Apply a linear martingale estimating function to find estimates of the
parameters of a process solution of a stochastic differential equation.
}
\usage{
linear.mart.ef(X, drift, sigma, a1, a2, guess, lower, upper, 
      c.mean, c.var) 
}
\arguments{
  \item{X}{a ts object containing a sample path of an sde.}
  \item{drift}{an expression for the drift coefficient; see details.}
  \item{sigma}{an expression for the diffusion coefficient; see details.}
  \item{a1, a2}{weights or instruments.}
  \item{c.mean}{expressions for the conditional mean.}
  \item{c.var}{expressions for the conditional variance.}
  \item{guess}{initial value of the parameters; see details.}
  \item{lower}{lower bounds for the parameters; see details.}
  \item{upper}{upper bounds for the parameters; see details.}
}
\details{
The function \code{linear.mart.ef} minimizes a linear martingale
estimating function that is a particular case of the polynomial
martingale estimating functions.
}
\value{
  \item{x}{a vector of estimates}
}
\references{
Bibby, B., Soerensen, M. (1995)  Martingale estimating functions for discretely 
observed diffusion processes,  \emph{Bernoulli}, 1, 17-39.
}
\author{Stefano Maria Iacus}
\examples{
set.seed(123)
d <- expression(-1 * x)
s <- expression(1) 
x0 <- rnorm(1,sd=sqrt(1/2))
sde.sim(X0=x0,drift=d, sigma=s,N=1000,delta=0.1) -> X
 
d <- expression(-theta * x)
  
linear.mart.ef(X, d, s, a1=expression(-x), lower=0, upper=Inf,
  c.mean=expression(x*exp(-theta*0.1)), 
  c.var=expression((1-exp(-2*theta*0.1))/(2*theta)))
}
\keyword{ts}
