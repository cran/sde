\name{gmm}
\alias{gmm}
\title{Generalized Method of Moments estimator}
\description{Implementation of the estimator of the Generalized Method of Moments by Hansen.
}
\usage{
gmm(X, u, dim, guess, lower, upper,maxiter=30,tol1=1e-3,tol2=1e-3)
}
\arguments{
  \item{X}{a ts object containg a sample path of a sde.}
  \item{u}{a function of `x', `y' and `theta'. See details.}
  \item{dim}{dimension of parameter space. See details.}
  \item{guess}{initial value of the parameters. See details.}
  \item{lower}{lower bounds for the parameters. See details.}
  \item{upper}{upper bounds for the parameters. See details.}
  \item{tol1}{tolerance for parameters. See details.}
  \item{tol2}{tolerance for Q1. See details.}
  \item{maxiter}{maximum number of iterations at the second stage. See details.}
}
\details{
The function \code{gmm} minimizes at the first stage
the  function \code{Q(theta) = t(Gn(theta)) * Gn(theta)} with respect to
\code{theta}, where \code{Gn(theta) =  mean(u(X[i+1], X[i], theta))}.
Then a  matrix of weights \code{W} is obtained by inverting an estimate
of the long-run covariance and the quadratic function 
\code{Q1(theta) = t(Gn(theta)) * W *  Gn(theta)} with starting value \code{theta1} (the solution
at the first stage).  The second stage is iterated until the first of these
conditions verifies: 1) number of iterations reaches \code{maxiter}; 2) euclidean
distance between theta1 and theta2 < tol1; 3) \code{Q1 < tol2}. 

The function \code{u} must be a function of \code{(u,y,theta)} and should
return a vector of the same length of the parameter space. The sanity checks
are left to the user.
}
\references{
Hansen, L.P. (1982) Large Sample Properties of Generalized Method of Moments Estimators, 
\emph{Econometrica}, 50(4),  1029-1054.
}
\value{
  \item{x}{a list with parameter estimates, the value of \code{Q1} at the minimum, the Hessian}
}
\author{Stefano Maria Iacus}
\note{This package is a companion to the book \emph{Simulation and Inference
for Stochastic Differential Equation}, Springer, NY.
}
\examples{
\dontrun{
alpha <- 0.5
beta <- 0.2
sigma <- sqrt(0.05)
true <- c(alpha, beta, sigma)
names(true) <- c("alpha", "beta", "sigma")

x0 <- rsCIR(1,theta=true)
set.seed(123)
sde.sim(X0=x0,model="CIR",theta=true,N=500000,delta=0.001) -> X
X <- window(X, deltat=0.1)
DELTA = deltat(X)
n <- length(X) 
X <- window(X, start=n*DELTA*0.5)
plot(X)

u <- function(x, y, theta){
  c.mean <- theta[1]/theta[2] + (y-theta[1]/theta[2])*exp(-theta[2]*DELTA)
  c.var <- ((y*theta[3]^2 * (exp(-theta[2]*DELTA)-exp(-2*theta[2]*DELTA))/theta[2] +  
				+theta[1]*theta[3]^2*(1-exp(-2*theta[2]*DELTA))/(2*theta[2]^2)))  
  cbind(x-c.mean,y*(x-c.mean), c.var-(x-c.mean)^2, y*(c.var-(x-c.mean)^2))  
}

CIR.lik <- function(theta1,theta2,theta3) {
 n <- length(X)
 dt <- deltat(X)
 -sum(dcCIR(x=X[2:n], Dt=dt, x0=X[1:(n-1)], theta=c(theta1,theta2,theta3), 
   log=TRUE))
}

fit <- mle(CIR.lik, start=list(theta1=.1,  theta2=.1,theta3=.3), 
    method="L-BFGS-B",lower=c(0.001,0.001,0.001), upper=c(1,1,1))
# maximum likelihood estimates
coef(fit)


gmm(X,u, guess=as.numeric(coef(fit)), lower=c(0,0,0), upper=c(1,1,1))

true
}}
\keyword{ts}
