\name{SIMloglik}
\alias{SIMloglik}
\title{Pedersen's approximation of the likelihood}
\description{Pedersen's approximation 
of the likelihood of a process solution of a stochastic differential equation.
This function is useful to calculate approximated maximum likelihood estimators
when the transition density of the process is not known. It is computationally
intensive.
}
\usage{
SIMloglik(X, theta, d, s,  M=10000, N=2, log=TRUE)
}
\arguments{
  \item{X}{a ts object containing a sample path of an sde.}
  \item{theta}{vector of parameters.}
  \item{d,s}{drift and diffusion coefficients; see details.}
  \item{log}{logical; if TRUE, the log-likelihood is returned.}
  \item{N}{number of subintervals; see details.}
  \item{M}{number of Monte Carlo simulations, which should be an even number; see details.}
}
\details{
The function \code{SIMloglik} returns the simulated log-likelihood obtained by 
Pedersen's method.
The functions \code{s} and \code{d} are the drift and diffusion
coefficients with arguments \code{(t,x,theta)}.
}
\value{
  \item{x}{a number}
}
\author{Stefano Maria Iacus}
\references{
Pedersen, A. R. (1995) A new approach to maximum likelihood estimation for 
stochastic differential equations based on discrete observations,
\emph{Scand. J. Statist.}, 22, 55-71.
}
\examples{
\dontrun{
set.seed(123)
d <- expression(-1*x)
s <- expression(2) 
sde.sim(drift=d, sigma=s,N=50,delta=0.01) -> X

S <- function(t, x, theta) sqrt(theta[2])
B <- function(t, x, theta) -theta[1]*x

true.loglik <- function(theta) {
 DELTA <- deltat(X)
 lik <- 0
 for(i in 2:length(X))
  lik <- lik + dnorm(X[i], mean=X[i-1]*exp(-theta[1]*DELTA), 
   sd = sqrt((1-exp(-2*theta[1]*DELTA))*
              theta[2]/(2*theta[1])),TRUE)
 lik  
}

xx <- seq(-10,10,length=20)
sapply(xx, function(x) true.loglik(c(x,4))) -> py
sapply(xx, function(x) EULERloglik(X,c(x,4),B,S)) -> pz
sapply(xx, function(x) SIMloglik(X,c(x,4),B,S,M=10000,N=5)) -> pw

plot(xx,py,type="l",xlab=expression(beta),
   ylab="log-likelihood",ylim=c(0,15)) # true
lines(xx,pz, lty=2) # Euler
lines(xx,pw, lty=3) # Simulated
}}
\keyword{ts}
