\name{sdeDiv}
\alias{sdeDiv}
\title{Phi-Divergences test for diffusion processes}
\description{Phi-Divergences test for diffusion processes.
}
\usage{
sdeDiv(X, theta1, theta0, phi= expression( -log(x) ), C.phi, K.phi, 
       b, s, b.x, s.x, s.xx, B, B.x, H, S, guess, \dots)
}
\arguments{
  \item{X}{a ts object containing a sample path of an sde.}
  \item{theta1}{a vector parameters for the hypothesis H1. If not given, \code{theta1} is estimated from the data.}
  \item{theta0}{a vector parameters for the hypothesis H0.}
  \item{phi}{an expression containing the phi function of the phi-divergence.}
  \item{C.phi}{the value of first derivtive of \code{phi} at point 1. If not given, it is calculated within this function.}
  \item{K.phi}{the value of second derivative of \code{phi} at point 1. If not given, it is calculated within this function.}
  \item{b}{drift coefficient of the model as a function of \code{x} and \code{theta}.}
  \item{s}{diffusion coefficient of the model as a function of \code{x} and \code{theta}.}
  \item{b.x}{partial derivative of \code{b} as a function of \code{x} and \code{theta}.}
  \item{s.x}{partial derivative of \code{s} as a function of \code{x} and \code{theta}.}
  \item{s.xx}{second-order partial derivative of \code{s} as a function of \code{x} and \code{theta}.}
  \item{B}{initial value of the parameters; see details.}
  \item{B.x}{partial derivative of \code{B} as a function of \code{x} and \code{theta}.}
  \item{H}{function of \code{(x,y)},  the integral of \code{B/s}; optional.}
  \item{S}{function of \code{(x,y)}, the integral of \code{1/s}; optional.}
  \item{guess}{initial value for the parameters to be estimated; optional.}
  \item{\dots}{passed to the \code{optim} function; optional.}
}
\details{
The \code{sdeDiv} estimate the phi-divergence for diffusion processes defined as
\code{D(theta1, theta0) = phi( f(theta1)/f(theta0) )} where \code{f} is the
likelihood function of the process. This function uses the Dacunha-Castelle 
and Florens-Zmirou approximation of the likelihood for \code{f}.

The parameter \code{theta1} is supposed to be the value of the true MLE estimator 
or the minimum contrast estimator of the parameters in the model. If missing
or \code{NULL} and \code{guess} is specified, \code{theta1} is estimated using the 
minimum contrast estimator derived from the locally Gaussian approximation
of the density. If both \code{theta1} and \code{guess} are missing, nothing can
be calculated.

The function always calculates the likelihood ratio test and the p-value of the
test statistics.
In some cases, the p-value of the phi-divergence test statistics is obtained by simulation. In such
a case, the \code{out$est.pval} is set to \code{TRUE}

Dy default \code{phi} is set to \code{-log(x)}. In this case the phi-divergence and the
likelihood ratio test are equivalent (e.g. phi-Div = LRT/2)

For more informations on phi-divergences for discretely observed diffusion processes
see the references.

If missing, \code{B} is calculated as \code{b/s - 0.5*s.x} provided that \code{s.x}
is not missing.

If missing, \code{B.x} is calculated as \code{b.x/s - b*s.x/(s^2)-0.5*s.xx}, provided
that \code{b.x}, \code{s.x}, and \code{s.xx} are not missing.

If missing, both \code{H} and \code{S} are evaluated numerically.
}

\references{
Dacunha-Castelle, D., Florens-Zmirou, D. (1986) Estimation of the coefficients 
of a diffusion from discrete observations,  \emph{Stochastics}, 19, 263-284.

De Gregorio, A., Iacus, S.M. (2008) Divergences Test Statistics for 
Discretely Observed Diffusion Processes, Journal of Statistical Planning and Inference, 140(7), 1744-1753, \doi{10.1016/j.jspi.2009.12.029}.
}
\value{
  \item{x}{a list containing the value of the divergence, its pvalue, the likelihood ratio test
  statistics and its p-value}
}
\author{Stefano Maria Iacus}
\examples{
\dontrun{
set.seed(123)
theta0 <- c(0.89218*0.09045,0.89218,sqrt(0.032742))
theta1 <- c(0.89218*0.09045/2,0.89218,sqrt(0.032742/2))

# we test the true model against two competing models
b <- function(x,theta) theta[1]-theta[2]*x
b.x <- function(x,theta)  -theta[2]

s <- function(x,theta) theta[3]*sqrt(x)
s.x <- function(x,theta) theta[3]/(2*sqrt(x))
s.xx <- function(x,theta) -theta[3]/(4*x^1.5)


X <- sde.sim(X0=rsCIR(1, theta1), N=1000, delta=1e-3, model="CIR", 
 theta=theta1)

sdeDiv(X=X, theta0 = theta0,   b=b, s=s, b.x=b.x, s.x=s.x, 
 s.xx=s.xx, method="L-BFGS-B", 
 lower=rep(1e-3,3), guess=c(1,1,1))

sdeDiv(X=X, theta0 = theta1,   b=b, s=s, b.x=b.x, s.x=s.x, 
 s.xx=s.xx, method="L-BFGS-B", 
 lower=rep(1e-3,3), guess=c(1,1,1))

lambda <- -1.75
myphi <-  expression( (x^(lambda+1) -x - lambda*(x-1))/(lambda*(lambda+1)) )

sdeDiv(X=X, theta0 = theta0,   phi = myphi, b=b, s=s, b.x=b.x, 
 s.x=s.x, s.xx=s.xx, method="L-BFGS-B", 
 lower=rep(1e-3,3), guess=c(1,1,1))

sdeDiv(X=X, theta0 = theta1,   phi = myphi, b=b, s=s, b.x=b.x, 
 s.x=s.x, s.xx=s.xx, method="L-BFGS-B", 
 lower=rep(1e-3,3), guess=c(1,1,1))
}
}
\keyword{ts}
