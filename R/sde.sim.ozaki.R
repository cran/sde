"sde.sim.ozaki" <-
function(X0, t0, Dt, N, d1, d1.x, s1){
   X <- numeric(N+1)
   B <- function(x) d1(1,x)
   Bx <- function(x) d1.x(1,x)
   S <- s1(1,1)   
   X[1] <- X0
   for(i in 2:(N+1)){
    x <- X[i-1]
    Kx <- log(1+B(x)*(exp(Bx(x)*Dt)-1)/(x*Bx(x)))/Dt
    Ex <- x + B(x)/Bx(x)*(exp(Bx(x)*Dt)-1)
    Vx <- S^2 * (exp(2*Kx*Dt) -1)/(2*Kx)
    X[i] <- rnorm(1, mean=Ex, sd=sqrt(Vx)) 
   }
 X
}
