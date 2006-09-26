\name{sde.sim}
\alias{sde.sim}
\title{Simulation of Stochastic Differential Equation}
\description{
Generic interface to different methods of simulations of solutions to stochastic
differential equations.
}
\usage{
sde.sim(t0 = 0, T = 1, X0 = 1, N = 100, delta, drift, sigma, drift.x, 
sigma.x, drift.xx, sigma.xx, drift.t, method = c("euler", "milstein", "KPS", 
"milstein2","cdist", "ozaki", "shoji"), alpha = 0.5, eta = 0.5, pred.corr = T,
rcdist=NULL, theta=NULL, model = c("CIR", "VAS", "OU", "BS"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t0}{time origin}
  \item{T}{horizon of simulation}
  \item{X0}{initial value of the process}
  \item{N}{number of simulation steps}
  \item{delta}{time-step of the simulation}
  \item{drift}{drift coeffcient: a expression of two variables \code{t} and \code{x}}
  \item{sigma}{diffusion coeffcient: a expression of two variables \code{t} and \code{x}}
  \item{drift.x}{partial derivative of drift coeffcient wrt to \code{x}: a function of two variables \code{t} and \code{x}}
  \item{sigma.x}{partial derivative of diffusio coefficient wrt to \code{x}: a function of two variables \code{t} and \code{x}}
  \item{drift.xx}{second partial derivative of drift coefficient wrt to \code{x}: a function of two variables \code{t} and \code{x}}
  \item{sigma.xx}{second partial derivative of diffusion coefficient wrt to \code{x}: a function of two variables \code{t} and \code{x}}
  \item{drift.t}{partial derivative of drift coefficient wrt to \code{t}: a function of two variables \code{t} and \code{x}}
  \item{method}{method of simulation, see details}
  \item{alpha}{weight \code{alpha} of the predictor-corrector scheme}
  \item{eta}{weight \code{eta} of the predictor-corrector scheme}
  \item{pred.corr}{boolean: wheter to apply the predictor-correct adjustment. See details.}
  \item{rcdist}{a function which is a random number generator from the conditional distribution of the process. See details.}
  \item{theta}{vector of parameters for \code{cdist}. See details.}
  \item{model}{model from which to simulate. See details.}
}
\details{
The function returns a \code{ts} object of length \code{N+1}, i.e. \code{X0} and
the new \code{N} simulated values. If \code{delta} is not specified, then \code{delta = (T-t0)/N}.
If \code{delta} is specified, then \code{N}  values of the solution of the sde are generated and
the time horizon \code{T} is adjusted to be \code{N * delta}.

If any of \code{drift.x}, \code{drift.xx}, \code{drift.t},
\code{sigma.x} and \code{sigma.xx} are not specified,
then numerical derivation is attempted when needed.

If \code{sigma} is not specified, it is assumed to the constant function \code{1}

The \code{method} of simulation can be one among: \code{euler}, \code{milstein},
\code{milstein2}, \code{KPS}, \code{cdist}, \code{ozaki} and \code{shoji}. 
No assumption on the coefficients or on \code{cdist} is checked: the user is 
responbile for using the right method for the process object of simulation.

The \code{model} is one among: \code{CIR}: Cox-Ingersoll-Ross, \code{VAS}: Vasicek,
\code{OU} Ornstein-Uhlenbeck, \code{BS}: Black and Scholes. 
No assumption on the coefficients \code{theta} is checked: the user is responbile 
for using the right ones.

If \code{method} is \code{cdist} then the process is simulated according to its
known conditional distribution. The random generator \code{rcdist} must be a 
function of \code{n}: the number of random numbers, \code{dt}: the time lag, 
\code{x}: the value of the process at time \code{t} - \code{dt} and the
vector of parameters \code{theta}.
}
\value{
  \item{x}{returns and invisible \code{ts} object}
}
\references{See Chapter 2 of the book}
\author{Stefano Maria Iacus}
\note{This package is a companion to the book `Simulation and Inference
for Stochastic Differential Equation, Springer, NY.
}
\examples{
# Ornstein-Uhlenbeck process
# beta = 5
# sigma = 3.5
set.seed(123)
d <- expression(-5 * x)
s <- expression(3.5) 
sde.sim(X0=10,drift=d, sigma=s) -> X
plot(X,main="Ornstein-Uhlenbeck")

# Cox-Ingersoll-Ross (CIR-1)
# a = 3
# theta = 2
# lambda = 5
set.seed(123)
d <- expression( 3*(2-x) ) 
s <- expression( 5*sqrt(x) ) 
sde.sim(X0=10,drift=d, sigma=s) -> X
plot(X,main="Cox-Ingersoll-Ross")

# Cox-Ingersoll-Ross using the conditional distribution
rCIR <- function(n=1, dt, x, theta){
   cc <- (2*theta[1])/((1-exp(-theta[1]*dt))*theta[3]^2)
   df <- 4*theta[1]*theta[2]/theta[3]^2
   ncp <- 2*cc*x*exp(-theta[1]*dt)
   rchisq(n, df=df, ncp=ncp)/(2*cc)
}

sde.sim(X0=10, theta=c(3, 2, 5), rcdist=rCIR, method="cdist") -> X
plot(X, main="Cox-Ingersoll-Ross")

sde.sim(X0=10, theta=c(3, 2, 5), model="CIR") -> X
plot(X, main="Cox-Ingersoll-Ross")
}
\keyword{datagen}
\keyword{ts}
