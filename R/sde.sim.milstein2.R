"sde.sim.milstein2" <-
function(X0,  t0, dt, N, 
          d1, d1.x, d1.xx, s1, s1.x, s1.xx){
   return( .Call("sde_sim_milstein2",  X0,  t0, dt, as.integer(N), 
              d1, d1.x, d1.xx, s1, s1.x, s1.xx, .GlobalEnv, PACKAGE="sde") )
}

