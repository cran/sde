#.First.lib <- function(lib, pkg) library.dynam("sde", pkg, lib) 

.noGenerics <- TRUE

.onLoad <- function(libname, pkgname)
{
# cat("Companion package to the book:\n")
# cat("Simulation and Inference for Stochastic Differential Equations With R Examples, Springer NY, (2008)\n")
 
 library.dynam("sde", pkgname, libname) 
}
