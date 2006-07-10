"sde.sim.euler" <-
function(X0, t0, dt, N, d1, s1, s1.x, alpha, eta, pred.corr){
 return( .Call("sde_sim_euler",  X0, t0, dt, as.integer(N), d1, s1, s1.x, 
              alpha, eta, as.logical(pred.corr), .GlobalEnv, PACKAGE="sde") ) 
}

