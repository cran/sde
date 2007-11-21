\name{HPloglik}
\alias{HPloglik}
\title{Ait-Sahalia  Hermite polynomial expansion  approximation 
of the likelihood}
\description{Ait-Sahalia Hermite polynomial expansion and Euler approximation 
of the likelihood of a process solution of a stochastic differential equation.
These functions are useful to calculate approximated maximum likelihood estimators
when the transition density of the process is not known.
}
\usage{
HPloglik(X, theta, M, F, s, log=TRUE)
}
\arguments{
  \item{X}{a ts object containing a sample path of an sde.}
  \item{theta}{vector of parameters.}
  \item{M}{list of derivatives; see details.}
  \item{F}{the transform function; see details.}
  \item{s}{drift and diffusion coefficient; see details.}
  \item{log}{logical; if TRUE, the log-likelihood is returned.}
}
\details{
The function \code{HPloglik} returns the Hermite polynomial approximation of the
likelihood of a diffusion process transformed to have a unitary diffusion
coefficient. The function \code{F} is the transform function, and
\code{s} is the original diffusion coefficient. The list of functions
\code{M} contains the transformed drift in \code{M[[1]]} and the
subsequent six derivatives in \code{x} of \code{M[[1]]}. The functions
\code{F}, \code{s}, and \code{M} have arguments \code{(t,x,theta)}.
}
\value{
  \item{x}{a number}
}
\references{ 
Ait-Sahalia, Y. (1996) Testing Continuous-Time Models of the Spot Interest Rate,  
\emph{Review of Financial Studies}, 9(2), 385-426.
}
\author{Stefano Maria Iacus}
\examples{
set.seed(123)
d <- expression(-1*x)
s <- expression(2) 
sde.sim(drift=d, sigma=s) -> X

M0 <- function(t, x, theta) -theta[1]*x
M1 <- function(t, x, theta) -theta[1]
M2 <- function(t, x, theta) 0
M3 <- function(t, x, theta) 0
M4 <- function(t, x, theta) 0
M5 <- function(t, x, theta) 0
M6 <- function(t, x, theta) 0
mu <- list(M0, M1, M2, M3, M4, M5, M6)

F <- function(t, x, theta) x/sqrt(theta[2])
S <- function(t, x, theta) sqrt(theta[2])

true.loglik <- function(theta) {
 DELTA <- deltat(X)
 lik <- 0
 for(i in 2:length(X))
  lik <- lik + dnorm(X[i], mean=X[i-1]*exp(-theta[1]*DELTA), 
   sd = sqrt((1-exp(-2*theta[1]*DELTA))*theta[2]/
             (2*theta[1])),TRUE)
 lik  
}

xx <- seq(-3,3,length=100)
sapply(xx, function(x) HPloglik(X,c(x,4),mu,F,S)) -> px
sapply(xx, function(x) true.loglik(c(x,4))) -> py

plot(xx,px,type="l",xlab=expression(beta),ylab="log-likelihood")
lines(xx,py, lty=3) # true
}
\keyword{ts}
