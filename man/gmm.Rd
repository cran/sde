\name{gmm}
\alias{gmm}
\title{Generalized method of moments estimator}
\description{Implementation of the estimator of the generalized method of moments by Hansen.
}
\usage{
gmm(X, u, dim, guess, lower, upper, maxiter=30, tol1=1e-3, 
   tol2=1e-3)
}
\arguments{
  \item{X}{a ts object containing a sample path of an sde.}
  \item{u}{a function of \code{x}, \code{y}, and \code{theta} and \code{DELTA}; see details.}
  \item{dim}{dimension of parameter space; see details.}
  \item{guess}{initial value of the parameters; see details.}
  \item{lower}{lower bounds for the parameters; see details.}
  \item{upper}{upper bounds for the parameters; see details.}
  \item{tol1}{tolerance for parameters; see details.}
  \item{tol2}{tolerance for Q1; see details.}
  \item{maxiter}{maximum number of iterations at the second stage; see details.}
}
\details{
The function \code{gmm} minimizes at the first stage
the  function \code{Q(theta) = t(Gn(theta)) * Gn(theta)} with respect to
\code{theta}, where \code{Gn(theta) =  mean(u(X[i+1], X[i], theta))}.
Then a  matrix of weights \code{W} is obtained by inverting an estimate
of the long-run covariance and the quadratic function 
\code{Q1(theta) = t(Gn(theta)) * W *  Gn(theta)} with starting value \code{theta1} (the solution
at the first stage).  The second stage is iterated until the first of these
conditions verifies: (1) that the number of iterations reaches \code{maxiter}; (2) that the Euclidean
distance between \code{theta1} and \code{theta2 < tol1}; (3) that \code{Q1 < tol2}. 

The function \code{u} must be a function of \code{(u,y,theta,DELTA)} and should
return a vector of the same length as the dimension of the parameter space. The sanity checks
are left to the user.
}
\references{
Hansen, L.P. (1982) Large Sample Properties of Generalized Method of Moments Estimators, 
\emph{Econometrica}, 50(4),  1029-1054.
}
\value{
  \item{x}{a list with parameter estimates, the value of \code{Q1} at the minimum, and the Hessian}
}
\author{Stefano Maria Iacus}
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

u <- function(x, y, theta, DELTA){
  c.mean <- theta[1]/theta[2] + 
             (y-theta[1]/theta[2])*exp(-theta[2]*DELTA)
  c.var <- ((y*theta[3]^2 * 
     (exp(-theta[2]*DELTA)-exp(-2*theta[2]*DELTA))/theta[2] +  
	  theta[1]*theta[3]^2*
	  (1-exp(-2*theta[2]*DELTA))/(2*theta[2]^2)))  
  cbind(x-c.mean,y*(x-c.mean), c.var-(x-c.mean)^2, 
        y*(c.var-(x-c.mean)^2))  
}

CIR.lik <- function(theta1,theta2,theta3) {
 n <- length(X)
 dt <- deltat(X)
 -sum(dcCIR(x=X[2:n], Dt=dt, x0=X[1:(n-1)], 
   theta=c(theta1,theta2,theta3), log=TRUE))
}

fit <- mle(CIR.lik, start=list(theta1=.1,  theta2=.1,theta3=.3), 
    method="L-BFGS-B",lower=c(0.001,0.001,0.001), upper=c(1,1,1))
# maximum likelihood estimates
coef(fit)


gmm(X,u, guess=as.numeric(coef(fit)), lower=c(0,0,0), 
    upper=c(1,1,1))

true
}}
\keyword{ts}
