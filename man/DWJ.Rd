\name{DWJ}
\alias{DWJ}
\docType{data}
\title{Weekly closings of the Dow-Jones industrial average}
\description{
This dataset contains the
weekly closings of the Dow-Jones industrial average in the period July 1971--August 1974.
These data were proposed to test
change-point estimators. There are 162 data, and the main evidence found by several authors 
is that a change in the variance occurred around the
third week of March 1973. 
}
\usage{data(DWJ)}
\references{
Hsu, D.A. (1977) Tests for variance shift at an unknown time point, \emph{Appl. Statist.}, 26(3), 279-284.

Hsu, D.A. (1979) Detecting shifts of parameter in gamma sequences with applications to stock price
and air traffic flow analysis, \emph{Journal American Stat. Ass.}, 74(365), 31-40.
}
\examples{
data(DWJ)
ret <- diff(DWJ)/DWJ[-length(DWJ)]

par(mfrow=c(2,1))
par(mar=c(3,3,2,1))
plot(DWJ,main="Dow-Jones closings",ylab="",type="p")
plot(ret,main="Dow-Jones returns",ylab="",type="p")

cp <- cpoint(ret)
cp
abline(v=cp$tau0,lty=3)

cp <- cpoint(window(ret,end=cp$tau0))
cp
abline(v=cp$tau0,lty=3)
}
\keyword{datasets}
