/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Code of this package: Copyright (C) 2006 S. M. Iacus
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 * Exports
 *	sde_sim_xxx(...)
 *
 * to be called as  .C(.)  in ../R/sde.sim.xxx.R
 * where xxx is one among "euler", "milstein", "milstei2", "KPS"
 */



#include <R.h>
#include <Rmath.h>
#include <R_ext/Boolean.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Complex.h>

#define max(a, b) (a > b ? a : b) 
#define min(a, b) (a < b ? a : b) 



/* Adapted from the code in man/doc/R-exts 
 * used to evaluate drift and diffusion coefficients
 */
double feval(double t, double x, SEXP f, SEXP rho);


/* interfaces to sde.sim.xxx(.R) */
/* parameters: do not apply to all functions

   x0: starting point
   t0: starting time
   delta: step of simulation
   N: number of simulated values (after x0)
   d: drift coefficient
   dx: partial derivative wrt to x of d
   dxx: second partial derivative wrt to x of d
   s: diffusion coefficient
   sx: partial derivative wrt to x of s
   sxx: second partial derivative wrt to x of s
   eta: weight in eta in the predictor-corrector adjustment
   alpha: weight in alpha in the predictor-corrector adjustment
   corr: apply predictor-correct adjustment?
   rho: the environtment on which to evaluate d, dx, dxx, s, sx, sxx
   Z: vector of pseudo-random gaussian numbers
   U: vector of pseudo-random numbers

 */ 
  
SEXP sde_sim_euler(SEXP x0, SEXP t0, SEXP delta, SEXP N, 
                   SEXP d, SEXP s, SEXP sx, 
				   SEXP eta, SEXP alpha, SEXP corr, SEXP rho);

SEXP sde_sim_milstein(SEXP x0, SEXP t0, SEXP delta, SEXP N,
                      SEXP d, SEXP s, SEXP sx, SEXP rho);

SEXP sde_sim_milstein2(SEXP x0, SEXP t0, SEXP delta, SEXP N,
                      SEXP d, SEXP dx, SEXP dxx,
					  SEXP s, SEXP sx, SEXP sxx, SEXP rho);

SEXP sde_sim_KPS(SEXP x0, SEXP t0, SEXP delta, SEXP N,
                      SEXP d, SEXP dx, SEXP dxx,
					  SEXP s, SEXP sx, SEXP sxx,
					  SEXP Z, SEXP U, SEXP rho);

/* Euler scheme */
SEXP sde_sim_euler(SEXP x0, SEXP t0, SEXP delta, SEXP N, 
                   SEXP d, SEXP s, SEXP sx, 
				   SEXP eta, SEXP alpha, SEXP corr, SEXP rho)
{
  SEXP X;
  double  Y1, Y2, T1, T2;
  double DELTA, ETA, ALPHA, *tm;
  double sdt, Z, tmp, d1, d2;
  Rboolean CORR;
  int i, n;
  
  if(!isNumeric(x0)) error("`x0' must be numeric");
  if(!isNumeric(t0)) error("`t0' must be numeric");
  if(!isNumeric(delta)) error("`delta' must be numeric");
  if(!isInteger(N)) error("`N' must be integer");

  if(!isFunction(d)) error("`d' must be a function");
  if(!isFunction(s)) error("`s' must be a function");
  if(!isFunction(sx)) error("`sx' must be a function");

  if(!isNumeric(eta)) error("`eta' must be numeric");
  if(!isNumeric(alpha)) error("`alpha' must be numeric");
  if(!isLogical(corr)) error("`corr' must be logical");

  if(!isEnvironment(rho)) error("`rho' must be an environment");
 
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(eta = AS_NUMERIC(eta));
  PROTECT(alpha = AS_NUMERIC(alpha));
  PROTECT(corr = AS_LOGICAL(corr));
  
  n = *INTEGER_POINTER(N);

  PROTECT(X = NEW_NUMERIC(n+1));
  REAL(X)[0] = *NUMERIC_POINTER(x0);

  T1 = *NUMERIC_POINTER(t0);
  DELTA = *NUMERIC_POINTER(delta);
  ETA = *NUMERIC_POINTER(eta);
  ALPHA = *NUMERIC_POINTER(alpha);
  CORR = *LOGICAL_POINTER(corr);

  sdt = sqrt(DELTA);

  Y1 = REAL(X)[0];
  
  GetRNGstate();
  if(CORR==TRUE){
   for(i=1; i<= n+1; i++){
    Z = rnorm(0,sdt);
    T2 = T1 + DELTA;
	tmp = REAL(X)[i-1];    
	Y2 = tmp + feval(T1,tmp,d, rho)*DELTA + feval(T1,tmp,s, rho)*Z;
	d1 = feval(T2,Y2,d,rho) - ETA*feval(T2,Y2, s, rho)*feval(T2,Y2,sx,rho);
	d2 = feval(T2,tmp,d,rho) - ETA*feval(T2,tmp, s, rho)*feval(T2,tmp,sx,rho);
	REAL(X)[i] = tmp + (ALPHA*d1 + (1.0-ALPHA)*d2)*DELTA +     
              (ETA * feval(T2,Y2,s,rho) + (1.0-ETA)*feval(T1,Y1,s,rho))*Z;
	
	Y1 = Y2;
	T1 = T2;  
   }   
  } else {
    for(i=1; i<= n+1; i++){
     Z = rnorm(0,sdt);  
	 T1 = T1 + DELTA; 
	 tmp = REAL(X)[i-1]; 
     REAL(X)[i] = tmp + feval(T1,tmp,d, rho)*DELTA + feval(T1,tmp,s, rho)*Z;
	}
  }

  PutRNGstate();

  UNPROTECT(7);
  return(X);
}


/* milstein scheme */
SEXP sde_sim_milstein(SEXP x0, SEXP t0, SEXP delta, SEXP N,
                      SEXP d, SEXP s, SEXP sx, SEXP rho)
{
  double T, DELTA;
  double sdt, Z, tmp, D, S, Sx;
  int i, n;
  SEXP X;
  
  if(!isNumeric(x0)) error("`x0' must be numeric");
  if(!isNumeric(t0)) error("`t0' must be numeric");
  if(!isNumeric(delta)) error("`delta' must be numeric");
  if(!isInteger(N)) error("`N' must be integer");  
  if(!isFunction(d)) error("`d' must be a function");
  if(!isFunction(s)) error("`s' must be a function");
  if(!isFunction(sx)) error("`sx' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
 
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(N = AS_INTEGER(N));
  
  T = *NUMERIC_POINTER(t0);
  n = *INTEGER_POINTER(N);
  DELTA = *NUMERIC_POINTER(delta);

  PROTECT(X = NEW_NUMERIC(n+1));
  REAL(X)[0] = *NUMERIC_POINTER(x0);
  sdt = sqrt(DELTA);

  GetRNGstate();
  for(i=1; i<= n+1; i++){
   Z = rnorm(0,sdt);
   T = T + DELTA;    
   tmp = REAL(X)[i-1];
   D = feval(T,tmp,d,rho);
   S = feval(T,tmp,s,rho);
   Sx = feval(T,tmp,sx,rho);
   REAL(X)[i] = tmp + D*DELTA + S*Z + 0.5*S*Sx*(Z*Z-DELTA);
  }   
  PutRNGstate();

  UNPROTECT(5);
  return(X);
}

/* second milstein scheme */
SEXP sde_sim_milstein2(SEXP x0, SEXP t0, SEXP delta, SEXP N,
                      SEXP d, SEXP dx, SEXP dxx,
					  SEXP s, SEXP sx, SEXP sxx, SEXP rho)
{
  double T, DELTA;
  double sdt, Z, tmp, D, Dx, Dxx, S, Sx, Sxx;
  int i, n;
  SEXP X;
  
  if(!isNumeric(x0)) error("`x0' must be numeric");
  if(!isNumeric(t0)) error("`t0' must be numeric");
  if(!isNumeric(delta)) error("`delta' must be numeric");
  if(!isInteger(N)) error("`N' must be integer");  
  if(!isFunction(d)) error("`d' must be a function");
  if(!isFunction(dx)) error("`dx' must be a function");
  if(!isFunction(dxx)) error("`dxx' must be a function");
  if(!isFunction(s)) error("`s' must be a function");
  if(!isFunction(sx)) error("`sx' must be a function");
  if(!isFunction(sxx)) error("`sxx' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
 
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(N = AS_INTEGER(N));
  
  T = *NUMERIC_POINTER(t0);
  n = *INTEGER_POINTER(N);
  DELTA = *NUMERIC_POINTER(delta);

  PROTECT(X = NEW_NUMERIC(n+1));
  REAL(X)[0] = *NUMERIC_POINTER(x0);
  sdt = sqrt(DELTA);

  GetRNGstate();
  for(i=1; i<= n+1; i++){
   Z = rnorm(0,sdt);
   T = T + DELTA;    
   tmp = REAL(X)[i-1];
   D = feval(T,tmp,d,rho);
   Dx = feval(T,tmp,dx,rho);
   Dxx = feval(T,tmp,dxx,rho);
   S = feval(T,tmp,s,rho);
   Sx = feval(T,tmp,sx,rho);
   Sxx = feval(T,tmp,sxx,rho);
   REAL(X)[i] = tmp + D*DELTA + S*Z + 0.5*S*Sx*(Z*Z-DELTA)+
			 pow(DELTA,1.5)*(0.5*D*Sx + 0.5*Dx*S + 0.25*S*Sxx)*Z +
			 DELTA*DELTA*(0.5*D*Dx+ 0.25*Dxx*S*S);
  }   
  PutRNGstate();

  UNPROTECT(5);
  return(X);
}

/* KPS scheme */
SEXP sde_sim_KPS(SEXP x0, SEXP t0, SEXP delta, SEXP N,
                      SEXP d, SEXP dx, SEXP dxx,
					  SEXP s, SEXP sx, SEXP sxx,
					  SEXP Z, SEXP U, SEXP rho)
{
  double T, DELTA, z, u;
  double sdt, tmp;
  double D, Dx, Dxx, S, Sx, Sxx;
  int i, n;
  SEXP X;
  
  if(!isNumeric(x0)) error("`x0' must be numeric");
  if(!isNumeric(t0)) error("`t0' must be numeric");
  if(!isNumeric(delta)) error("`delta' must be numeric");
  if(!isInteger(N)) error("`N' must be integer");  
  
  if(!isFunction(d)) error("`d' must be a function");
  if(!isFunction(dx)) error("`dx' must be a function");
  if(!isFunction(dxx)) error("`dxx' must be a function");
  
  if(!isFunction(s)) error("`s' must be a function");
  if(!isFunction(sx)) error("`sx' must be a function");
  if(!isFunction(sxx)) error("`sxx' must be a function");
  
  if(!isNumeric(Z)) error("`Z' must be numeric");
  if(!isNumeric(U)) error("`U' must be numeric");

  if(!isEnvironment(rho)) error("`rho' must be an environment");
 
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(Z = AS_NUMERIC(Z));
  PROTECT(U = AS_NUMERIC(U));
  PROTECT(N = AS_INTEGER(N));
  
  T = *NUMERIC_POINTER(t0);
  n = *INTEGER_POINTER(N);
  DELTA = *NUMERIC_POINTER(delta);

  PROTECT(X = NEW_NUMERIC(n+1));
  REAL(X)[0] = *NUMERIC_POINTER(x0);
  sdt = sqrt(DELTA);

  for(i=1; i<= n+1; i++){
    T = T + DELTA;    
    tmp = REAL(X)[i-1];
   
	D = feval(T,tmp,d,rho);
    Dx = feval(T,tmp,dx,rho);
    Dxx = feval(T,tmp,dxx,rho);
    S = feval(T,tmp,s,rho);
    Sx = feval(T,tmp,sx,rho);
    Sxx = feval(T,tmp,sxx,rho);
    z = REAL(Z)[i-1];
	u = REAL(U)[i-1];

     REAL(X)[i] = tmp + D * DELTA + S * z +
             0.5 * S * Sx * (z*z-DELTA) + 
			 S * Dx * u +
			 0.5 * (D * Dx + 0.5 * S*S * Dxx) * DELTA*DELTA +
			 (D * Sx + 0.5 * S*S * Sxx) * (z * DELTA - u) +
			 0.5 * S * (Sx*Sx + S*Sxx) * (z*z/3.0 - DELTA) * z;
  }   

  UNPROTECT(7);
  return(X);
}


static R_CMethodDef R_CDef[] = {
   {"sde_sim_euler", (DL_FUNC)&sde_sim_euler, 11},
   {"sde_sim_milstein", (DL_FUNC)&sde_sim_milstein, 8},
   {"sde_sim_milstein2", (DL_FUNC)&sde_sim_milstein2, 11},
   {"sde_sim_KPS", (DL_FUNC)&sde_sim_KPS, 13},
   {NULL, NULL, 0},
};

void
R_init_ifs(DllInfo *info)
{
    R_registerRoutines(info, R_CDef, NULL, NULL, NULL);
}




/* feval:
   ------------
   t   : time variable
   x   : space variable
   f   : a SEXP to a R function
   rho : the environment `f' is going to be evaluated
   
   on return:
   ----------
   
   the value of f(t,x)   
*/


double feval(double t, double x, SEXP f, SEXP rho)
{
    double val= 0.0;
    SEXP R_fcall, tpar, xpar; 
	
	PROTECT(tpar = allocVector(REALSXP, 1));
	PROTECT(xpar = allocVector(REALSXP, 1));
    REAL(tpar)[0] = t;
    REAL(xpar)[0] = x;
   
	PROTECT(R_fcall = allocList(3));
	SETCAR(R_fcall, f);
	SET_TYPEOF(R_fcall, LANGSXP);

	SETCADR(R_fcall, tpar);
	SETCADDR(R_fcall, xpar);
    val = *NUMERIC_POINTER(eval(R_fcall, rho));
    UNPROTECT(3);
		
    return(val);

}

