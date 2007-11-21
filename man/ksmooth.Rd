\name{ksmooth}
\alias{ksdrift}
\alias{ksdiff}
\alias{ksdens}
\title{Nonparametric invariant density, drift, and diffusion coefficient estimation}
\description{Implementation of simple Nadaraya-Watson nonparametric 
estimation of drift and diffusion
coefficient, and plain kernel density estimation of the invariant density
for a  one-dimensional diffusion process.}
\usage{
ksdrift(x, bw, n = 512)
ksdiff(x, bw, n = 512)
ksdens(x, bw, n = 512)
}
\arguments{
  \item{x}{a \code{ts} object.}
  \item{bw}{bandwidth.}
  \item{n}{number of points in which to calculate the estimates.}
}
\details{
These functions return the nonparametric estimate of the drift or
diffusion coefficients for data \code{x} using the Nadaraya-Watson estimator 
for diffusion processes.


\code{ksdens} returns the density estimates of the invariant density.


If not provided, the bandwidth \code{bw}
is calculated using Scott's rule (i.e., 
\code{bw = len^(-1/5)*sd(x)}) where \code{len=length(x)} 
is the number of observed points of the diffusion path.
}
\value{
  \item{val}{an invisible list of \code{x} and \code{y} coordinates 
  and an object of class \code{density} in the case of invariant
  density estimation}
}
\author{Stefano Maria Iacus}
\references{
Ait-Sahalia, Y. (1996) Nonparametric pricing of interest rate 
derivative securities, \emph{Econometrica}, 64, 527-560.

Bandi, F., Phillips, P. (2003) Fully nonparametric estimation of  
scalar diffusion models, \emph{Econometrica}, 71, 241-283.

Florens-Zmirou, D. (1993) On estimating the diffusion coefficient 
from discrete observations,  \emph{Journal of Applied Probability}, 30, 790-804.
}
\examples{
set.seed(123)
theta <- c(6,2,1)
X <- sde.sim(X0 = rsCIR(1, theta), model="CIR", theta=theta,
     N=1000,delta=0.1)

b <- function(x)
 theta[1]-theta[2]*x

sigma <- function(x)
 theta[3]*sqrt(x)
  
minX <- min(X)
maxX <- max(X)

par(mfrow=c(3,1))
curve(b,minX,maxX)
lines(ksdrift(X),lty=3)

curve(sigma,minX, maxX)
lines(ksdiff(X),lty=3)

f <-function(x) dsCIR(x, theta)
curve(f,minX,maxX)
lines(ksdens(X),lty=3)
}
\keyword{datagen}
\keyword{ts}
