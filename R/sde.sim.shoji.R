"sde.sim.shoji" <-
function(X0, t0, dt, N, d1, d1.x, d1.xx, d1.t, s1){
   X <- numeric(N+1)
   S <- s1(1,1)   
   X[1] <- X0
   for(i in 2:(N+1)){
    x <- X[i-1]
    Lx <- d1.x(dt,x)
	Mx <- S^2 * d1.xx(dt,x)/2 + d1.t(dt,x)
	Ex <- (x + d1(dt,x)*(exp(Lx*dt)-1)/Lx + 
	       Mx*(exp(Lx*dt) -1 -Lx*dt)/Lx^2) 
    Vx <- S^2*(exp(2*Lx*dt)-1)/(2*Lx)
    X[i] <- rnorm(1, mean=Ex, sd=sqrt(Vx)) 
   }
 X
}

