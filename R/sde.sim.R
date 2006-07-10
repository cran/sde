"sde.sim" <-
function(t0=0, T=1, X0=1, N=100, delta,
         drift, sigma, drift.x, sigma.x, drift.xx, sigma.xx,
         scheme=c("euler","milstein","KPS","milstein2"),
         alpha=0.5, eta=0.5, pred.corr=T){

  if(missing(drift))
   stop("please specify al least the drift coefficient of the SDE")

  if(missing(sigma))
   sigma <- expression(1)

  if(!is.expression(drift) | !is.expression(sigma))
   stop("coefficients must be expressions in `t' an `x'") 

  scheme <- match.arg(scheme)
  
  if(pred.corr==F){
   alpha <- 0
   eta <- 0
   sigma.x <- NULL
  }

  needs.sx <- FALSE;
  needs.dx <- FALSE;
  needs.sxx <- FALSE;
  needs.dxx <- FALSE;

  
  if(scheme=="milstein") needs.sx = TRUE
  if((scheme=="euler" & pred.corr==T)) 
   needs.sx = TRUE
  if(scheme == "KPS" | scheme == "milstein2") {
   needs.sx <- TRUE
   needs.dx <- TRUE
   needs.sxx <- TRUE
   needs.dxx <- TRUE
  }
  
  
  if(needs.sx & missing(sigma.x)){
    cat("sigma.x not provided, attempting symbolic derivation.\n")
    sigma.x <- D(sigma,"x")
   }
  
  if(needs.dx & missing(drift.x)){
	cat("drift.x not provided, attempting symbolic derivation.\n")
	drift.x <- D(drift,"x")
   }

   if(needs.dxx & missing(drift.xx)){
	cat("drift.xx not provided, attempting symbolic derivation.\n")
	drift.xx <- D(D(drift,"x"),"x")
   }
   

   if(needs.sxx & missing(sigma.xx)){
	cat("sigma.xx not provided, attempting symbolic derivation.\n")
	sigma.xx <- D(D(sigma,"x"),"x")
   }

 
  d1 <- function(t,x)  eval(drift)
  d1.x <- function(t,x) eval(drift.x)
  d1.xx <- function(t,x) eval(drift.x)
  s1 <- function(t,x) eval(sigma)
  s1.x <- function(t,x) eval(sigma.x)
  s1.xx <- function(t,x) eval(sigma.xx)

  if(t0<0 | T<0)
   stop("please use positives times!")
  
  if(missing(delta)){
    t <- seq(t0,T, length=N+1)
  } else {
   t <- c(t0,t0+cumsum(rep(delta,N)))
   T <- t[N+1]
   warning("T set to =",T,"\n") 
  }

   dt <- (T-t0)/N

  if(scheme=="euler")
   X <- sde.sim.euler(X0, t0, dt, N, d1, s1, s1.x, alpha, eta, pred.corr)  
  
  if(scheme=="milstein")
   X <- sde.sim.milstein(X0,  t0, dt, N, d1, s1, s1.x)
  
  if(scheme=="milstein2")
   X <- sde.sim.milstein2(X0,  t0, dt, N, d1, d1.x, d1.xx, s1, s1.x, s1.xx)

  if(scheme=="KPS"){
  	require(MASS)
	Sigma <- matrix(c(dt, 0.5*dt^2, 0.5*dt^2, 1/3*dt^3),2,2)
 	tmp <- mvrnorm(N, c(0,0), Sigma)
	Z <- tmp[,1]
	U <- tmp[,2]
    X <- sde.sim.KPS(X0,  t0, dt, N, d1, d1.x, d1.xx, 
	           s1, s1.x, s1.xx, Z, U)
  }

  X <- ts(X, start=t0, deltat=dt)
  invisible(X)
}

