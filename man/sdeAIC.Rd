\name{sdeAIC}
\alias{sdeAIC}
\title{Akaike's Infromation Criterion for diffusion processes}
\description{Implementation of the AIC statistics for diffuzion processes.
}
\usage{
sdeAIC(X, theta, b, s, b.x, s.x, s.xx, B, B.x, H, S)
}
\arguments{
  \item{X}{a ts object containg a sample path of a sde.}
  \item{theta}{a vector or estimates of the parameters.}
  \item{b}{drift coefficient of the model as a function of \code{x} and \code{theta}}
  \item{s}{diffusion coefficient of the model as a function of \code{x} and \code{theta}}
  \item{b.x}{partial derivative of \code{b} as a function of \code{x} and \code{theta}}
  \item{s.x}{partial derivative of \code{s} as a function of \code{x} and \code{theta}}
  \item{s.xx}{second order partial derivative of \code{s} as a function of \code{x} and \code{theta}}
  \item{B}{initial value of the parameters. See details.}
  \item{B.x}{partial deriivative of \code{B} as a function of \code{x} and \code{theta}}
  \item{H}{function of \code{(x,y)}. It is the integral of \code{B/s}. Optional.}
  \item{S}{function of \code{(x,y)}. It is the integral of \code{1/s}. Optional.}
}
\details{
the \code{sdeAIC} evaluates the AIC statistics for diffusion processes using
Dacunha-Castelle and Florens-Zmirou approximation of the likelihood.

If missing, \code{B} is calculated as \code{b/s - 0.5*s.x} provided that \code{s.x}
is not missing.

If missing, \code{B.x} is calculated as \code{b.x/s - b*s.x/(s^2)-0.5*s.xx} provided
that \code{b.x}, \code{s.x} and \code{s.xx} are not missing.

If missing, both \code{H} and \code{S} are evaluated numerically.
}

\references{
Dacunha-Castelle, D., Florens-Zmirou, D. (1986) Estimation of the coefficients 
of a diffusion from discrete observations,  \emph{Stochastics}, 19, 263-284.

Uchida, M., Yoshida, N.  (2005) AIC for ergodic diffusion processes from 
discrete observations, preprint MHF 2005-12, march 2005, \emph{Faculty of 
Mathematics, Kyushu University, Fukuoka, Japan}.
}
\value{
  \item{x}{the value of the AIC statistics}
}
\author{Stefano Maria Iacus}
\note{This package is a companion to the book \emph{Simulation and Inference
for Stochastic Differential Equation}, Springer, NY.
}
\examples{
set.seed(123)
dri <- expression(-(x-10))
dif <- expression(2) 
sde.sim(X0=10,drift=dri, sigma=dif,N=1000,delta=0.1) -> X

b <- function(x,theta) -theta[1]*(x-theta[2])
b.x <- function(x,theta)  -theta[1]+0*x

# the true model
s <- function(x,theta) theta[3]*sqrt(x)
s.x <- function(x,theta) theta[3]/(2*sqrt(x))
s.xx <- function(x,theta) -theta[3]/(4*x^1.5)
sdeAIC(X, c(1,10,2), b, s, b.x, s.x, s.xx)

# not the true model
s <- function(x,theta) sqrt(theta[3]*+theta[4]*x)
s.x <- function(x,theta) theta[4]/(2*sqrt(theta[3]+theta[4]*x))
s.xx <- function(x,theta) -theta[4]^2/(4*(theta[3]+theta[4]*x)^1.5)
sdeAIC(X, c(1,10,2,4), b, s, b.x, s.x, s.xx)
}
\keyword{ts}
