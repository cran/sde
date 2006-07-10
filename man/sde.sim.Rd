\name{sde.sim}
\alias{sde.sim}
\title{Simulation of Stochastic Differential Equation}
\description{
Generic interface to different scheme of simulations of solutions to stochastic
differential equations.
}
\usage{
sde.sim(t0 = 0, T = 1, X0 = 1, N = 100, delta, drift, sigma, drift.x, 
sigma.x, drift.xx, sigma.xx, scheme = c("euler", "milstein", "KPS", "milstein2"), 
alpha = 0.5, eta = 0.5, pred.corr = T)
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
  \item{sigma.x}{partial derivative of diffusio coeffcient wrt to \code{x}: a function of two variables \code{t} and \code{x}}
  \item{drift.xx}{second partial derivative of drift coeffcient wrt to \code{x}: a function of two variables \code{t} and \code{x}}
  \item{sigma.xx}{second partial derivative of diffusion coeffcient wrt to \code{x}: a function of two variables \code{t} and \code{x}}
  \item{scheme}{scheme of simulation, see details}
  \item{alpha}{weight \code{alpha} of the predictor-corrector scheme}
  \item{eta}{weight \code{eta} of the predictor-corrector scheme}
  \item{pred.corr}{boolean: wheter to apply the predictor-correct adjustment. See details.}
}
\details{
The function returns a \code{ts} object of length \code{N+1}, i.e. \code{X0} and
the new \code{N} simulated values. If \code{delta} is not specified, then \code{delta = (T-t0)/N}.
If \code{delta} is specified, then \code{N}  values of the solution of the sde are generated and
the time horizon \code{T} is adjusted to be \code{N * delta}.

If \code{d.x} and/or \code{d.xx} and/or \code{s.x} and/or \code{s.xx} are not specified,
then numerical derivation is attempted.

If \code{s} is not specified, it is assumed to the constant function \code{1}

The \code{scheme} of simulation can be one among: \code{euler}, \code{milstein},
\code{milstein2} or \code{KPS}. No assumption on the coefficients is checked: the
user is responbile for using the right scheme for the process object of simulation.


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
}
\keyword{datagen}
\keyword{ts}
