"sde.sim.KPS" <-
function(X0,  t0, dt, N, d1, d1.x, d1.xx, s1, s1.x, s1.xx, Z, U){
    return( .Call("sde_sim_KPS",  X0,  t0, dt, as.integer(N), d1, d1.x, d1.xx, 
	           s1, s1.x, s1.xx, Z, U, .GlobalEnv, PACKAGE="sde") )
}

