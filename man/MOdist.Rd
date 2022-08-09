\name{MOdist}
\alias{MOdist}
\title{Markov Operator distance for clustering diffusion processes.}
\description{Markov Operator distance for clustering diffusion processes.}
\usage{
MOdist(x, M=50, rangeval=range(x, na.rm=TRUE, finite = TRUE))
}
\arguments{
  \item{x}{one or multi-dimensional time series.}
  \item{M}{number of splines bases used to approximate the Markov Operator.}
  \item{rangeval}{a vector containing lower and upper limit. Default is the range of \code{x}.}
}
\details{
This function return a lower triangular dist object to be further used 
in cluster analysis (see examples below).

If \code{x} is a one-dimensional time series, the output is the scalar 0, not a \code{dist} object.

If \code{x} has less than 2 observations, \code{NA} is returned.

If time series \code{x} contains missing data, then \code{x} is converted to a \code{zoo} object and
missing data are imputed by interpolation.
}
\value{
  \item{X}{a \code{dist} object}
}
\author{Stefano Maria Iacus}
\references{
De Gregorio, A. Iacus, S.M. (2008) Clustering of discretely observed diffusion processes, Computational Statistics and Data Analysis, 54(12), 598-606, \doi{10.1016/j.csda.2009.10.005}.
}
\examples{
\dontrun{
data(quotes)

plot(quotes)

d <- MOdist(quotes)
cl <- hclust( d )
groups <- cutree(cl, k=4)

cmd <- cmdscale(d)
plot( cmd, col=groups)
text( cmd, labels(d) , col=groups)

plot(quotes, col=groups)

plot(quotes, col=groups,ylim=range(quotes))
}
}
\keyword{datagen}
\keyword{ts}
