"sde.sim.ea" <-
function(t0=0, T=1, X0=1, N=100, delta, drift,
 drift.x, k1, k2, phi, max.psi = 1000, rh, A){
   if(missing(drift.x)){ 
     cat("drift.x not provided, attempting symbolic derivation.\n")
     drift.x <- D(drift,"x")
   }

   if(t0<0 | T<0)
    stop("please use positives times!")
 
   d1 <- function(x) eval(drift)
   d1.x <- function(x)  eval(drift.x)
   psi <- function(x) 0.5*d1(x)^2 + 0.5*d1.x(x)

   if(missing(k1)){
    cat("k1 missing, trying numerical minimization...") 
    k1 <- optimize(psi, c(0, max.psi))$obj
    cat(sprintf("(k1=%5.3f)\n",k1))
   }
   if(missing(k2)){
    cat("k2 missing, trying numerical maximization...") 
    k2 <- optimize(psi, c(0, max.psi),max=TRUE)$obj
    cat(sprintf("(k2=%5.3f)\n",k2))
  }
   
  if(missing(phi))
    phi <- function(x) 0.5*d1(x) + 0.5*d1.x(x) - k1
  else
   phi <- function(x) eval(phi)

  M <- k2-k1
  if(M==0)
   stop("`k1' = `k2' probably due to numerical maximization")
   
  if(missing(delta)){
    t <- seq(t0, T, length=N+1)
  } else {
   t <- c(t0,t0+cumsum(rep(delta,N)))
   T <- t[N+1]
   cat(sprintf("T set to %4.3f\n",T)) 
  }

   dt <- (T-t0)/N
   if(dt>1/M)
    stop(sprintf("discretization step greater than 1/(k2_k1)"))
   
   if(missing(A))
     A <- function(x) integrate(d1, 0, x)

  if(missing(rh)){
   rh <- function(){
    h <- function(x) exp(A(x) - x^2/(2*dt))
    f <- function(x) h(x)/dnorm(x,sqrt(dt))
    maxF <- optimize(f,c(-3*dt, 3*dt),max=TRUE)$obj
    while(1){
     y <- rnorm(1)
     if( runif(1) < f(y)/maxF )
      return(y)
   }
  }
 }

  x0 <- X0
  X <- numeric(N)
  X[1] <- X0
  rej <- 0
  j <- 1
  while(j <= N){
   y <- x0+rh()
   k <- rpois(1,M*dt)
   if(k>0){
    t <- runif(k)*dt
    v <- runif(k)*M
    idx <- order(t)
	t <- c(0, t[idx], dt)
    v <- v[idx]

	Dt <- t[2:(k+2)] - t[1:(k+1)]
	W <- c(0,cumsum(sqrt(Dt) * rnorm(k+1))) 
    Y <- x0 + W -(W[k+2] -y+x0)*t/dt
    if( prod(phi(Y[2:(k+1)]) <= v) == 1){
     j <- j+1
     x0 <- Y[k+2]
     X[j] <-  Y[k+2]
    } else {
     rej <- rej +1
   }
 }
}
cat(sprintf("rejection rate: %5.3f\n",rej/(N+rej)))
return(invisible(ts(X,start=t0, delta=dt)))
}

