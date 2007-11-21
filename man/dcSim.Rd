\name{dcSim}
\alias{dcSim}
\title{Pedersen's simulated transition density}
\description{
Simulated transition density \eqn{X(t) | X(t_0) = x_0}X(t) | X(t0) = x0 of a diffusion process based on
Pedersen's method.}
\usage{
dcSim(x0, x, t, d, s, theta, M=10000, N=10, log=FALSE)
}
\arguments{
  \item{x0}{the value of the  process at time \code{0}.}
  \item{x}{value in which to evaluate the conditional density.}
  \item{t}{lag or time.}
  \item{theta}{parameter of the process; see details.}
  \item{log}{logical; if TRUE, probabilities \eqn{p}{p} are given as \eqn{\log(p)}{log(p)}.}
  \item{d}{drift coefficient as a function; see details.}
  \item{s}{diffusion coefficient as a function; see details.}
  \item{N}{number of subintervals; see details.}
  \item{M}{number of Monte Carlo simulations, which should be an even number; see details.}
}
\details{
This function returns the value of the conditional density of
\eqn{X(t) | X(0) = x_0}{X(t) | X(0) = x0} at point \code{x}. 

The functions \code{d} and \code{s}, must be functions of \code{t}, 
\code{x}, and \code{theta}.
}
\value{
  \item{x}{a numeric vector}
}
\references{
Pedersen, A. R. (1995) A new approach to maximum likelihood estimation for 
stochastic differential equations based on discrete observations,
\emph{Scand. J. Statist.}, 22, 55-71.
}
\author{Stefano Maria Iacus}
\examples{
\dontrun{
d1 <- function(t,x,theta) theta[1]*(theta[2]-x)
s1 <- function(t,x,theta) theta[3]*sqrt(x)

from <- 0.08
x <- seq(0,0.2, length=100)
sle10 <- NULL
sle2 <- NULL
sle5 <- NULL
true <- NULL
set.seed(123)
for(to in x){
 sle2 <- c(sle2, dcSim(from, to, 0.5, d1, s1, 
    theta=c(2,0.02,0.15), M=50000,N=2))
 sle5 <- c(sle5, dcSim(from, to, 0.5, d1, s1, 
    theta=c(2,0.02,0.15), M=50000,N=5))
 sle10 <- c(sle10, dcSim(from, to, 0.5, d1, s1, 
    theta=c(2,0.02,0.15), M=50000,N=10))
 true <- c(true, dcCIR(to, 0.5, from, c(2*0.02,2,0.15)))
}

par(mar=c(5,5,1,1))
plot(x, true, type="l", ylab="conditional density")
lines(x, sle2, lty=4)
lines(x, sle5, lty=2)
lines(x, sle10, lty=3)
legend(0.15,20, legend=c("exact","N=2", "N=5", "N=10"), 
   lty=c(1,2,4,3))
}}
\keyword{datagen}
\keyword{ts}
