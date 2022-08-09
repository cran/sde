\name{quotes}
\alias{quotes}
\docType{data}
\title{Daily closings of 20 financial time series from 2006-01-03 to 2007-12-31}
\description{
This dataset contains the
daily closings of 20 assets from NYSE/NASDAQ. Quotes from 2006-01-03 to 2007-12-31.

It is an object of class \code{zoo}. Original data contained missing data and interpolated.
Used as example data to test the Markov Operator
distance for discretely observed diffusion processes.
}
\usage{data(quotes)}
\references{
De Gregorio, A. Iacus, S.M. (2008) Clustering of discretely observed diffusion processes, Computational Statistics and Data Analysis, 54(12), 598-606, \doi{10.1016/j.csda.2009.10.005}.
}
\examples{
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
\keyword{datasets}
