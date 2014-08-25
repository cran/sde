\name{sdeAIC}
\alias{sdeAIC}
\title{Akaike's information criterion for diffusion processes}
\description{Implementation of the AIC statistics for diffusion processes.
}
\usage{
sdeAIC(X, theta, b, s, b.x, s.x, s.xx, B, B.x, H, S, guess, 
       ...)
}
\arguments{
  \item{X}{a ts object containing a sample path of an sde.}
  \item{theta}{a vector or estimates of the parameters.}
  \item{b}{drift coefficient of the model as a function of \code{x} and \code{theta}.}
  \item{s}{diffusion coefficient of the model as a function of \code{x} and \code{theta}.}
  \item{b.x}{partial derivative of \code{b} as a function of \code{x} and \code{theta}.}
  \item{s.x}{partial derivative of \code{s} as a function of \code{x} and \code{theta}.}
  \item{s.xx}{second-order partial derivative of \code{s} as a function of \code{x} and \code{theta}.}
  \item{B}{initial value of the parameters; see details.}
  \item{B.x}{partial derivative of \code{B} as a function of \code{x} and \code{theta}.}
  \item{H}{function of \code{(x,y)},  the integral of \code{B/s}; optional.}
  \item{S}{function of \code{(x,y)}, the integral of \code{1/s}; optional.}
  \item{guess}{initial value for the parameters to be estimated; optional.}
  \item{...}{passed to the \code{optim} function; optional.}
}
\details{
The \code{sdeAIC} evaluates the AIC statistics for diffusion processes using
Dacunha-Castelle and Florens-Zmirou approximations of the likelihood.

The parameter \code{theta} is supposed to be the value of the true MLE estimator 
or the minimum contrast estimator of the parameters in the model. If missing
or \code{NULL} and \code{guess} is specified, \code{theta} is estimated using the 
minimum contrast estimator derived from the locally Gaussian approximation
of the density. If both \code{theta} and \code{guess} are missing, nothing can
be calculated.

If missing, \code{B} is calculated as \code{b/s - 0.5*s.x} provided that \code{s.x}
is not missing.

If missing, \code{B.x} is calculated as \code{b.x/s - b*s.x/(s^2)-0.5*s.xx}, provided
that \code{b.x}, \code{s.x}, and \code{s.xx} are not missing.

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
\examples{
\dontrun{
set.seed(123)
# true model generating data
dri <- expression(-(x-10))
dif <- expression(2*sqrt(x)) 
sde.sim(X0=10,drift=dri, sigma=dif,N=1000,delta=0.1) -> X

# we test the true model against two competing models
b <- function(x,theta) -theta[1]*(x-theta[2])
b.x <- function(x,theta)  -theta[1]+0*x

s <- function(x,theta) theta[3]*sqrt(x)
s.x <- function(x,theta) theta[3]/(2*sqrt(x))
s.xx <- function(x,theta) -theta[3]/(4*x^1.5)
# AIC for the true model
sdeAIC(X, NULL, b, s, b.x, s.x, s.xx, guess=c(1,1,1),
       lower=rep(1e-3,3), method="L-BFGS-B")

s <- function(x,theta) sqrt(theta[3]*+theta[4]*x)
s.x <- function(x,theta) theta[4]/(2*sqrt(theta[3]+theta[4]*x))
s.xx <- function(x,theta) -theta[4]^2/(4*(theta[3]+theta[4]*x)^1.5)
# AIC for competing model 1
sdeAIC(X, NULL, b, s, b.x, s.x, s.xx, guess=c(1,1,1,1),
       lower=rep(1e-3,4), method="L-BFGS-B")

s <- function(x,theta) (theta[3]+theta[4]*x)^theta[5]
s.x <- function(x,theta) 
           theta[4]*theta[5]*(theta[3]+theta[4]*x)^(-1+theta[5])
s.xx <- function(x,theta) (theta[4]^2*theta[5]*(theta[5]-1)
                 *(theta[3]+theta[4]*x)^(-2+theta[5]))
# AIC for competing model 2
sdeAIC(X, NULL, b, s, b.x, s.x, s.xx, guess=c(1,1,1,1,1),
       lower=rep(1e-3,5), method="L-BFGS-B")
}
}
\keyword{ts}
