\name{BM}
\alias{BM}
\alias{BBridge}
\alias{GBM}
\title{Brownian motion, Brownian Bridge and Geometric Brownian motion simulators}
\description{Brownian motion, Brownian Bridge and Geometric Brownian motion simulators}
\usage{
BBridge(x=0, y=0, t0=0, T=1, N=100)
BM(x=0, t0=0, T=1, N=100)
GBM(x=1, r=0, sigma=1, T=1, N=100)
}
\arguments{
  \item{x}{intial value of the process at time \code{t0}}
  \item{y}{terminal value of the process at time \code{T}}
  \item{t0}{initial time.}
  \item{r}{the interest rate of the GBM.}
  \item{sigma}{the volatility of the GBM.}
  \item{T}{final time.}
  \item{N}{number of intervals in which to split \code{[t0,T]}.}
}
\details{
These functions return an invisible \code{ts} object containing
a trajectory of the process calculated on grid of \code{N+1}
equidistant points between \code{t0} and \code{T}, i.e. 
\code{t[i] = t0 + (T-t0)*i/N}, \code{i in 0:N}.	\code{t0=0} for the
Geometric Brownian Motion.

The function \code{BBridge} returns a trajectory of the Brownian Bridge 
starting in \code{x} at time \code{t0} and
ending at \code{y} at time \code{T}, i.e. 
\code{(B(t), t0 <= t <= T | B(t_0)=x, B(T)=y)}

The function \code{BM} returns 
a trajectory of the translated 
Brownian Motion \code{(B(t), t>= t0 | B(t0)=x)}, 
i.e. \code{x+B(t-t0)}, for \code{t >= t0}. 
The standard Brownian motion is obtained
choosing \code{x=0} and \code{t0=0} (the default values).

The function \code{GBM} returns a trajectory of the Geometric Brownian Motion
starting at \code{x} at time \code{t0=0}, i.e. the process
\code{S(t) = x * exp((r-sigma^2/2)*t + sigma*B(t))}.
}
\value{
  \item{X}{an invisible \code{ts} object}
}
\author{Stefano Maria Iacus}
\note{This package is a companion to the book \emph{Simulation and Inference
for Stochastic Differential Equation}, Springer, NY.
}
\examples{
plot(BM())
plot(BBridge())
plot(GBM())
}
\keyword{datagen}
\keyword{ts}
