"sde.sim.cdist" <-
function(X0, t0, Dt, N, rcdist=NULL, theta=NULL){
   X <- numeric(N+1)
   X[1] <- X0
   for(i in 2:(N+1)){
    X[i] <- rcdist(1, Dt, X[i-1], theta)
   }
 X
}

