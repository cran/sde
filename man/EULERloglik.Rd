\name{EULERloglik}
\alias{EULERloglik}
\title{Euler approximation of the likelihood}
\description{Euler approximation 
of the likelihood of a process solution of a stochastic differential equation.
These functions are useful to calculate approximated maximum likelihood estimators
when the transition density of the process is not known.
}
\usage{
EULERloglik(X, theta, d, s, log = TRUE) 
}
\arguments{
  \item{X}{a ts object containing a sample path of an sde.}
  \item{theta}{vector of parameters.}
  \item{d,s}{drift and diffusion coefficients; see details.}
  \item{log}{logical; if TRUE, the log-likelihood is returned.}
}
\details{
The function \code{EULERloglik} returns the Euler approximation of the
log-likelihood. The functions \code{s} and \code{d} are the drift and diffusion
coefficients with arguments \code{(t,x,theta)}.
}
\value{
  \item{x}{a number}
}
\author{Stefano Maria Iacus}
\examples{
set.seed(123)
d <- expression(-1*x)
s <- expression(2) 
sde.sim(drift=d, sigma=s) -> X

S <- function(t, x, theta) sqrt(theta[2])
B <- function(t, x, theta) -theta[1]*x

true.loglik <- function(theta){
 DELTA <- deltat(X)
 lik <- 0
 for(i in 2:length(X))
  lik <- lik + dnorm(X[i], mean=X[i-1]*exp(-theta[1]*DELTA), 
  sd = sqrt((1-exp(-2*theta[1]*DELTA))*
            theta[2]/(2*theta[1])),TRUE)
 lik  
}

xx <- seq(-3,3,length=100)
sapply(xx, function(x) true.loglik(c(x,4))) -> py
sapply(xx, function(x) EULERloglik(X,c(x,4),B,S)) -> pz

# true likelihood
plot(xx,py,type="l",xlab=expression(beta),ylab="log-likelihood") 
lines(xx,pz, lty=2) # Euler
}
\keyword{ts}
