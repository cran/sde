\name{sde.sim}
\alias{sde.sim}
\title{Simulation of stochastic differential equation}
\description{
Generic interface to different methods of simulation of solutions to stochastic
differential equations.
}
\usage{
sde.sim(t0 = 0, T = 1, X0 = 1, N = 100, delta, drift, sigma, 
   drift.x, sigma.x, drift.xx, sigma.xx, drift.t, 
   method = c("euler", "milstein", "KPS", "milstein2", 
   "cdist","ozaki","shoji","EA"), 
   alpha = 0.5, eta = 0.5, pred.corr = T, rcdist = NULL, 
   theta = NULL, model = c("CIR", "VAS", "OU", "BS"),
   k1, k2, phi, max.psi = 1000, rh, A, M=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t0}{time origin.}
  \item{T}{horizon of simulation.}
  \item{X0}{initial value of the process.}
  \item{N}{number of simulation steps.}
  \item{M}{number of trajectories.}
  \item{delta}{time step of the simulation.}
  \item{drift}{drift coefficient: an expression of two variables \code{t} and \code{x}.}
  \item{sigma}{diffusion coefficient: an expression of two variables \code{t} and \code{x}.}
  \item{drift.x}{partial derivative of the drift coefficient w.r.t. \code{x}: a function of two variables \code{t} and \code{x}.}
  \item{sigma.x}{partial derivative of the diffusion coefficient w.r.t.  \code{x}: a function of two variables \code{t} and \code{x}.}
  \item{drift.xx}{second partial derivative of the drift coefficient w.r.t.  \code{x}: a function of two variables \code{t} and \code{x}.}
  \item{sigma.xx}{second partial derivative of the diffusion coefficient w.r.t. \code{x}: a function of two variables \code{t} and \code{x}.}
  \item{drift.t}{partial derivative of the drift coefficient w.r.t.  \code{t}: a function of two variables \code{t} and \code{x}.}
  \item{method}{method of simulation; see details.}
  \item{alpha}{weight \code{alpha} of the predictor-corrector scheme.}
  \item{eta}{weight \code{eta} of the predictor-corrector scheme.}
  \item{pred.corr}{boolean: whether to apply the predictor-correct adjustment; see details.}
  \item{rcdist}{a function that is a random number generator from the conditional distribution of the process; see details.}
  \item{theta}{vector of parameters for \code{cdist}; see details.}
  \item{model}{model from which to simulate; see details.}
  \item{k1}{lower bound for \code{psi(x)}; see details.}
  \item{k2}{upper bound for \code{psi(x)}; see details.}
  \item{phi}{the function \code{psi(x) - k1}.}
  \item{max.psi}{upper value of the support of \code{psi} to search for its maximum.}
  \item{rh}{the rejection function; see details.}
  \item{A}{\code{A(x)} is the integral of the \code{drift} between \code{0} and \code{x}.}
}
\details{
The function returns a \code{ts} object of length \code{N+1}; i.e., \code{X0} and
the new \code{N} simulated values if \code{M=1}.
For \code{M>1}, an \code{mts} (multidimensional \code{ts} object) is returned, which
means that \code{M} independent trajectories are simulated. 
If the initial value \code{X0} is not of the length \code{M}, the values are recycled
in order to have an initial vector of the correct length. 
If \code{delta} is not specified, then \code{delta = (T-t0)/N}.
If \code{delta} is specified, then \code{N}  values of the solution of the sde are generated and
the time horizon \code{T} is adjusted to be \code{N * delta}.

The function \code{psi} is \code{psi(x) = 0.5*drift(x)^2 + 0.5*drift.x(x)}.

If any of \code{drift.x}, \code{drift.xx}, \code{drift.t},
\code{sigma.x}, and \code{sigma.xx} are not specified,
then numerical derivation is attempted when needed.

If \code{sigma} is not specified, it is assumed to be the constant function \code{1}.

The \code{method} of simulation can be one among: \code{euler}, \code{KPS}, \code{milstein},
\code{milstein2},  \code{cdist}, \code{EA}, \code{ozaki}, and \code{shoji}. 
No assumption on the coefficients or on \code{cdist} is checked: the user is 
responsible for using the right method for the process object of simulation.

The \code{model} is one among: \code{CIR}: Cox-Ingersoll-Ross, \code{VAS}: Vasicek,
\code{OU} Ornstein-Uhlenbeck, \code{BS}: Black and Scholes. 
No assumption on the coefficient \code{theta} is checked: the user is responsible 
for using the right ones.

If the \code{method} is \code{cdist}, then the process is simulated according to its
known conditional distribution. The random generator \code{rcdist} must be a 
function of \code{n}, the number of random numbers; \code{dt}, the time lag; 
\code{x}, the value of the process at time \code{t} - \code{dt}; and the
vector of parameters \code{theta}.

For the exact algorithm method \code{EA}: if missing \code{k1} and \code{k2} as well 
as \code{A}, \code{rh} and \code{phi} are calculated numerically by the function.
}
\value{
  \item{x}{returns an invisible \code{ts} object}
}
\references{See Chapter 2 of the text.}
\author{Stefano Maria Iacus}
\examples{
# Ornstein-Uhlenbeck process
set.seed(123)
d <- expression(-5 * x)
s <- expression(3.5) 
sde.sim(X0=10,drift=d, sigma=s) -> X
plot(X,main="Ornstein-Uhlenbeck")

# Multiple trajectories of the O-U process
set.seed(123)
sde.sim(X0=10,drift=d, sigma=s, M=3) -> X
plot(X,main="Multiple trajectories of O-U")

# Cox-Ingersoll-Ross process
# dXt = (6-3*Xt)*dt + 2*sqrt(Xt)*dWt
set.seed(123)
d <- expression( 6-3*x ) 
s <- expression( 2*sqrt(x) ) 
sde.sim(X0=10,drift=d, sigma=s) -> X
plot(X,main="Cox-Ingersoll-Ross")

# Cox-Ingersoll-Ross using the conditional distribution "rcCIR"

set.seed(123)
sde.sim(X0=10, theta=c(6, 3, 2), rcdist=rcCIR, 
        method="cdist") -> X
plot(X, main="Cox-Ingersoll-Ross")

set.seed(123)
sde.sim(X0=10, theta=c(6, 3, 2), model="CIR") -> X
plot(X, main="Cox-Ingersoll-Ross")

# Exact simulation
set.seed(123)
d <- expression(sin(x))
d.x <- expression(cos(x)) 
A <- function(x) 1-cos(x)
sde.sim(method="EA", delta=1/20, X0=0, N=500, 
        drift=d, drift.x = d.x, A=A) -> X
plot(X, main="Periodic drift")
}
\keyword{datagen}
\keyword{ts}
