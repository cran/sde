\name{DBridge}
\alias{DBridge}
\title{Simulation of diffusion bridge}
\description{Simulation of diffusion bridge.}
\usage{
DBridge(x=0, y=0, t0=0, T=1, delta, drift, sigma, ...)
}
\arguments{
  \item{x}{initial value of the process at time \code{t0}.}
  \item{y}{terminal value of the process at time \code{T}.}
  \item{t0}{initial time.}
  \item{delta}{time step of the simulation.}
  \item{drift}{drift coefficient: an expression of two variables \code{t} and \code{x}.}
  \item{sigma}{diffusion coefficient: an expression of two variables \code{t} and \code{x}.}
  \item{T}{final time.}
  \item{...}{passed to the \code{sde.sim} function.}
}
\details{
The function returns a trajectory of the diffusion bridge 
starting at \code{x} at time \code{t0} and
ending at \code{y} at time \code{T}.


The function uses the \code{sde.sim} function to simulate the paths internally.
Refer to the \code{sde.sim} documentation for further information about the
argument ``\code{...}''
}
\value{
  \item{X}{an invisible \code{ts} object}
}
\author{Stefano Maria Iacus}
\references{
Bladt, M., Soerensen, M. (2007) Simple simulation of diffusion bridges 
with application to likelihood inference for diffusions, mimeo.
}
\seealso{
  \code{\link{sde.sim}}, \code{\link{BBridge}}
}
\examples{
d <- expression((3-x))
s <- expression(1.2*sqrt(x))
par(mar=c(3,3,1,1))
par(mfrow=c(2,1)) 
set.seed(123)
X <- DBridge(x=1.7,y=0.5, delta=0.01, drift=d, sigma=s)
plot(X)
X <- DBridge(x=1,y=5, delta=0.01, drift=d, sigma=s)
plot(X)
}
\keyword{datagen}
\keyword{ts}
