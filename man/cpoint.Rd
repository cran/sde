\name{cpoint}
\alias{cpoint}
\title{Volatility change-point estimator for diffusion processes}
\description{Volatility change-point estimator for diffusion processes based on least squares.}
\usage{
cpoint(x, mu, sigma)
}
\arguments{
  \item{x}{a \code{ts} object.}
  \item{mu}{a function of \code{x} describing the drift coefficient.}
  \item{sigma}{a function of \code{x} describing the diffusion coefficient.}
}
\details{
The function returns a list of  elements containing the discrete \code{k0} and continuous \code{tau0}
change-point instant, the estimated volatilities before (\code{theta1}) and after (\code{theta2}) the time change.
The model is assumed to be of the  form
\deqn{{\rm d}X_t = b(X_t){\rm d}t + \theta \sigma(X_t){\rm d}W_t}{dXt = b(Xt)dt + theta*sigma(Xt)dWt}
where \code{theta} = \code{theta1} for \code{t<=tau0} and \code{theta} = \code{theta2} otherwise.

If the drift coefficient is unknown, the  model
\deqn{{\rm d}X_t = b(X_t){\rm d}t + \theta {\rm d}W_t}{dXt = b(Xt)dt + \theta*dWt}
is considered and \code{b} is estimated nonparametrically.
}
\value{
  \item{X}{a list}
}
\author{Stefano Maria Iacus}
\examples{
tau0 <- 0.6
k0 <- ceiling(1000*tau0)
set.seed(123)
X1 <- sde.sim(X0=1, N=2*k0, t0=0, T=tau0, model="CIR", 
              theta=c(6,2,1))
X2 <- sde.sim(X0=X1[2*k0+1], N=2*(1000-k0), t0=tau0, 
   T=1, model="CIR", theta=c(6,2,3))

Y <- ts(c(X1,X2[-1]), start=0, deltat=deltat(X1))
X <- window(Y,deltat=0.01) 
DELTA <- deltat(X)
n <- length(X)

mu <- function(x) 6-2*x
sigma <- function(x) sqrt(x)

cp <- cpoint(X,mu,sigma)
cp
plot(X)
abline(v=tau0,lty=3)
abline(v=cp$tau0,col="red")

# nonparametric estimation
cpoint(X)
}
\keyword{datagen}
\keyword{ts}
