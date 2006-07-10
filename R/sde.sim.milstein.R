"sde.sim.milstein" <-
function(X0,  t0, dt, N, d1, s1, s1.x){
   return( .Call("sde_sim_milstein",  X0,  t0, dt, as.integer(N), d1, 
			  s1, s1.x, .GlobalEnv, PACKAGE="sde") )
}

