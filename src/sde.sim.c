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
double ftheta(double t, double x, SEXP theta, SEXP f, SEXP rho);

/* transition densities */



/* cHP:  approximated transition density. Internal function only.
		 to be called by HPloglik
 */
double cHP(double Delta, double mu0, double mu1, double mu2, double mu3,
						   double mu4, double mu5, double mu6, 
						   double z, double s);


/* HPloglik: likelihood approximation by Hermite poynomials 
   Ait-Sahalia, Y. (2002) "Maximum likelihood estimation of
   discretely sampled diffusions: a clsed form approximation 
   approach", Econometrica, 70, 223-262
*/
SEXP HPloglik(SEXP delta, SEXP X, SEXP theta, SEXP M0, SEXP M1,
                  SEXP M2, SEXP M3, SEXP M4, SEXP M5, SEXP M6, 
				  SEXP F, SEXP S, SEXP rho);

/* Euler log-likelihood */
SEXP EULERloglik(SEXP delta, SEXP X, SEXP theta, SEXP d, SEXP s, SEXP rho);


/* Pedersen's simulated transition density */
SEXP dcSim(SEXP X, SEXP Y, SEXP delta, SEXP d, SEXP s, SEXP theta, SEXP N, SEXP M, SEXP rho);

/* Pedersen's simulated likelihood */
SEXP SIMloglik(SEXP X, SEXP delta, SEXP d, SEXP s, SEXP theta, SEXP N, SEXP M, SEXP rho);

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
  
SEXP sde_sim_euler(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                   SEXP d, SEXP s, SEXP sx, 
				   SEXP eta, SEXP alpha, SEXP corr, SEXP rho);


SEXP sde_sim_milstein(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                      SEXP d, SEXP s, SEXP sx, SEXP rho);

SEXP sde_sim_milstein2(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                      SEXP d, SEXP dx, SEXP dxx,
					  SEXP s, SEXP sx, SEXP sxx, SEXP rho);

SEXP sde_sim_KPS(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                      SEXP d, SEXP dx, SEXP dxx,
					  SEXP s, SEXP sx, SEXP sxx,
					  SEXP Z, SEXP U, SEXP rho);

SEXP sde_sim_cdist(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                      SEXP cdist, SEXP rho);

SEXP sde_sim_ozaki(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                      SEXP d, SEXP dx, SEXP s, SEXP rho);

SEXP sde_sim_shoji(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                      SEXP d, SEXP dx, SEXP dxx, SEXP dt, SEXP s, SEXP rho);


SEXP sde_sim_euler(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                   SEXP d, SEXP s, SEXP sx, 
				   SEXP eta, SEXP alpha, SEXP corr, SEXP rho)
{
  SEXP X, Y1, Y2;
  double  *rY1, *rY2, T1, T2, *rX;
  double DELTA, ETA, ALPHA;
  double sdt, Z, tmp, d1, d2, *rx0;
  Rboolean CORR;
  int i, n, j, m;
  
  if(!isNumeric(x0)) error("`x0' must be numeric");
  if(!isNumeric(t0)) error("`t0' must be numeric");
  if(!isNumeric(delta)) error("`delta' must be numeric");
  if(!isInteger(N)) error("`N' must be integer");
  if(!isInteger(M)) error("`M' must be integer");

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
  
  n = *INTEGER(N);
  m = *INTEGER(M);

  PROTECT(Y1 = NEW_NUMERIC(m));
  PROTECT(Y2 = NEW_NUMERIC(m));
  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1));
  rX = REAL(X);
  rY1 = REAL(Y1);
  rY2 = REAL(Y2);
  rx0 = REAL(x0);  
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 
  T1 = *REAL(t0);
  DELTA = *REAL(delta);
  ETA = *REAL(eta);
  ALPHA = *REAL(alpha);
  CORR = *LOGICAL(corr);

  sdt = sqrt(DELTA);

  for(j=0; j<m; j++)
   rY1[j] = rX[j*(n+1)];
  
  GetRNGstate();
  if(CORR==TRUE){
   for(i=1; i< n+1; i++){
    T2 = T1 + DELTA;
    for(j=0; j<m; j++){
 	 Z = rnorm(0,sdt);
	 tmp = rX[i-1 +j*(n+1)];    
	 rY2[j] = tmp + feval(T1,tmp,d, rho)*DELTA + feval(T1,tmp,s, rho)*Z;
	 d1 = feval(T2,rY2[j],d,rho) - ETA*feval(T2,rY2[j], s, rho)*feval(T2,rY2[j],sx,rho);
	 d2 = feval(T2,tmp,d,rho) - ETA*feval(T2,tmp, s, rho)*feval(T2,tmp,sx,rho);
	 rX[i +j*(n+1)] = tmp + (ALPHA*d1 + (1.0-ALPHA)*d2)*DELTA +     
               (ETA * feval(T2,rY2[j],s,rho) + (1.0-ETA)*feval(T1,rY1[j],s,rho))*Z;	
 	 rY1[j] = rY2[j];
	}
	T1 = T2;  
   }   
  } else {
    for(i=1; i< n+1; i++){
	 T1 = T1 + DELTA; 
     for(j=0; j<m; j++){
	  Z = rnorm(0, sdt);  
	  tmp = rX[i + j*(n+1) - 1]; 
      rX[i + (n+1)*j] = tmp + feval(T1,tmp,d, rho)*DELTA + feval(T1,tmp,s, rho)*Z;
	  }
	}
  }

  PutRNGstate();

  UNPROTECT(9);
  return(X);
}


/* milstein method */
SEXP sde_sim_milstein(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                      SEXP d, SEXP s, SEXP sx, SEXP rho)
{
  double T, DELTA, *rX, *rx0;
  double sdt, Z, tmp, D, S, Sx;
  int i, n, j, m;
  SEXP X;
  
  if(!isNumeric(x0)) error("`x0' must be numeric");
  if(!isNumeric(t0)) error("`t0' must be numeric");
  if(!isNumeric(delta)) error("`delta' must be numeric");
  if(!isInteger(N)) error("`N' must be integer");  
  if(!isInteger(M)) error("`M' must be integer");  
  if(!isFunction(d)) error("`d' must be a function");
  if(!isFunction(s)) error("`s' must be a function");
  if(!isFunction(sx)) error("`sx' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
 
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(N = AS_INTEGER(N));
  
  T = *REAL(t0);
  n = *INTEGER(N);
  m = *INTEGER(M);

  DELTA = *REAL(delta);

  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1));
  rx0 = REAL(x0);  
  rX = REAL(X);
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 
  sdt = sqrt(DELTA);

  GetRNGstate();
  for(i=1; i< n+1; i++){
   T = T + DELTA;   
   for(j=0; j<m; j++){ 
    Z = rnorm(0,sdt);
	tmp = rX[i + j*(n+1) - 1]; 
    D = feval(T,tmp,d,rho);
    S = feval(T,tmp,s,rho);
    Sx = feval(T,tmp,sx,rho);
	rX[i + (n+1)*j] = tmp + D*DELTA + S*Z + 0.5*S*Sx*(Z*Z-DELTA);
   }
  }   
  PutRNGstate();

  UNPROTECT(5);
  return(X);
}

/* second milstein method */
SEXP sde_sim_milstein2(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                      SEXP d, SEXP dx, SEXP dxx,
					  SEXP s, SEXP sx, SEXP sxx, SEXP rho)
{
  double T, DELTA, *rX, *rx0;
  double sdt, Z, tmp, D, Dx, Dxx, S, Sx, Sxx;
  int i, n, j, m;
  SEXP X;
  
  if(!isNumeric(x0)) error("`x0' must be numeric");
  if(!isNumeric(t0)) error("`t0' must be numeric");
  if(!isNumeric(delta)) error("`delta' must be numeric");
  if(!isInteger(N)) error("`N' must be integer");  
  if(!isInteger(M)) error("`M' must be integer");  
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
  
  T = *REAL(t0);
  n = *INTEGER(N);
  m = *INTEGER(M);
  DELTA = *REAL(delta);

  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1));
  rx0 = REAL(x0);  
  rX = REAL(X);
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 

  sdt = sqrt(DELTA);

  GetRNGstate();
  for(i=1; i< n+1; i++){
   T = T + DELTA;
   for(j=0; j<m;j++){    
    Z = rnorm(0,sdt);
	tmp = rX[i + j*(n+1) - 1];  
    D = feval(T,tmp,d,rho);
    Dx = feval(T,tmp,dx,rho);
    Dxx = feval(T,tmp,dxx,rho);
    S = feval(T,tmp,s,rho);
    Sx = feval(T,tmp,sx,rho);
    Sxx = feval(T,tmp,sxx,rho);
    rX[i + j*(n+1)] = tmp + D*DELTA + S*Z + 0.5*S*Sx*(Z*Z-DELTA)+
			 pow(DELTA,1.5)*(0.5*D*Sx + 0.5*Dx*S + 0.25*S*Sxx)*Z +
			 DELTA*DELTA*(0.5*D*Dx+ 0.25*Dxx*S*S);
	}
  }   
  PutRNGstate();

  UNPROTECT(5);
  return(X);
}

/* KPS method */
SEXP sde_sim_KPS(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                      SEXP d, SEXP dx, SEXP dxx,
					  SEXP s, SEXP sx, SEXP sxx,
					  SEXP Z, SEXP U, SEXP rho)
{
  double T, DELTA, z, u, *rX, *rx0, *rZ, *rU;
  //  double sdt, tmp;
  double  tmp;
  double D, Dx, Dxx, S, Sx, Sxx;
  int i, n, j, m;
  SEXP X;
  
  if(!isNumeric(x0)) error("`x0' must be numeric");
  if(!isNumeric(t0)) error("`t0' must be numeric");
  if(!isNumeric(delta)) error("`delta' must be numeric");
  if(!isInteger(N)) error("`N' must be integer");  
  if(!isInteger(M)) error("`M' must be integer");  

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
  
  T = *REAL(t0);
  n = *INTEGER(N);
  m = *INTEGER(M);

  DELTA = *REAL(delta);

  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1));
  rx0 = REAL(x0);  
  rX = REAL(X);
  rZ = REAL(Z);
  rU = REAL(U);
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 

  //sdt = sqrt(DELTA);

  for(i=1; i< n+1; i++){
    T = T + DELTA;    
    for(j=0;j<m;j++){
     tmp = rX[i + j*(n+1) -1];
	 D = feval(T,tmp,d,rho);
     Dx = feval(T,tmp,dx,rho);
     Dxx = feval(T,tmp,dxx,rho);
     S = feval(T,tmp,s,rho);
     Sx = feval(T,tmp,sx,rho);
     Sxx = feval(T,tmp,sxx,rho);
     z = rZ[j*n + i-1];
	 u = rU[j*n + i-1];

     rX[i + j*(n+1)] = tmp + D * DELTA + S * z +
             0.5 * S * Sx * (z*z-DELTA) + 
			 S * Dx * u +
			 0.5 * (D * Dx + 0.5 * S*S * Dxx) * DELTA*DELTA +
			 (D * Sx + 0.5 * S*S * Sxx) * (z * DELTA - u) +
			 0.5 * S * (Sx*Sx + S*Sxx) * (z*z/3.0 - DELTA) * z;
	}
  }   

  UNPROTECT(7);
  return(X);
}

/* cdist method */
SEXP sde_sim_cdist(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                      SEXP cdist, SEXP rho)
{
  //  double T, DELTA, *rX, *rx0;
  double  DELTA, *rX, *rx0;
  int i, n, j, m;
  SEXP X;
  
  if(!isNumeric(x0)) error("`x0' must be numeric");
  if(!isNumeric(t0)) error("`t0' must be numeric");
  if(!isNumeric(delta)) error("`delta' must be numeric");
  if(!isInteger(N)) error("`N' must be integer");  
  if(!isInteger(M)) error("`M' must be integer");  
  if(!isFunction(cdist)) error("`cdist' must be a function");
  if(!isEnvironment(rho)) error("`rho' must be an environment");
 
  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(N = AS_INTEGER(N));
  
  //T = *REAL(t0);
  n = *INTEGER(N);
  m = *INTEGER(M);

  DELTA = *REAL(delta);

  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1));
  rx0 = REAL(x0);  
  rX = REAL(X);
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 

  GetRNGstate();
  for(i=1; i< n+1; i++)
   for(j=0; j<m; j++) 
	rX[i + (n+1)*j] = feval(DELTA, rX[i-1+j*(n+1)], cdist, rho);
   
  PutRNGstate();

  UNPROTECT(5);
  return(X);
}

/* Ozaki method */
SEXP sde_sim_ozaki(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                      SEXP d, SEXP dx, SEXP s, SEXP rho)
{
//    double T, DELTA, *rX, *rx0;
  double DELTA, *rX, *rx0;
  double tmp, D, Dx, S, Ex, Vx, Kx;
  int i, n, j, m;
  SEXP X;
  
  if(!isNumeric(x0)) error("`x0' must be numeric");
  if(!isNumeric(t0)) error("`t0' must be numeric");
  if(!isNumeric(delta)) error("`delta' must be numeric");
  if(!isInteger(N)) error("`N' must be integer");  
  if(!isInteger(M)) error("`M' must be integer");  
  if(!isFunction(d)) error("`d' must be a function");
  if(!isFunction(dx)) error("`dx' must be a function");
  if(!isNumeric(s)) error("`s' must be numeric");
  if(!isEnvironment(rho)) error("`rho' must be an environment");

  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(s = AS_NUMERIC(s));
  PROTECT(N = AS_INTEGER(N));
  
 // T = *REAL(t0);
  S = *REAL(s);
  n = *INTEGER(N);
  m = *INTEGER(M);
  DELTA = *REAL(delta);

  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1));
  rx0 = REAL(x0);  
  rX = REAL(X);
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 

  GetRNGstate();
  for(i=1; i< n+1; i++){
   for(j=0; j<m;j++){    
	tmp = rX[i + j*(n+1) - 1];  
    D = feval(1.0,tmp,d,rho);
    Dx = feval(1.0,tmp,dx,rho);
    Kx = log(1.0+D*(exp(Dx*DELTA)-1.0)/(tmp*Dx))/DELTA;
    Ex = tmp + D/Dx*(exp(Dx*DELTA)-1.0);
    Vx = S*sqrt((exp(2.0*Kx*DELTA) -1.0)/(2.0*Kx));    
	rX[i + j*(n+1)] = rnorm(Ex,Vx);
	}
  }   
  PutRNGstate();

  UNPROTECT(6);
  return(X);
}

/* Shoji method */
SEXP sde_sim_shoji(SEXP x0, SEXP t0, SEXP delta, SEXP N, SEXP M,
                      SEXP d, SEXP dx, SEXP dxx, SEXP dt, SEXP s, SEXP rho)
{
  //double T, DELTA, *rX, *rx0;
  double DELTA, *rX, *rx0;
  double tmp, D, Dx, Dxx, Dt, S, Ex, Vx, Mx;
  int i, n, j, m;
  SEXP X;
  
  if(!isNumeric(x0)) error("`x0' must be numeric");
  if(!isNumeric(t0)) error("`t0' must be numeric");
  if(!isNumeric(delta)) error("`delta' must be numeric");
  if(!isInteger(N)) error("`N' must be integer");  
  if(!isInteger(M)) error("`M' must be integer");  
  if(!isFunction(d)) error("`d' must be a function");
  if(!isFunction(dx)) error("`dx' must be a function");
  if(!isFunction(dxx)) error("`dxx' must be a function");
  if(!isFunction(dt)) error("`dt' must be a function");
  if(!isNumeric(s)) error("`s' must be numeric");
  if(!isEnvironment(rho)) error("`rho' must be an environment");

  PROTECT(x0 = AS_NUMERIC(x0));
  PROTECT(delta = AS_NUMERIC(delta));
  PROTECT(t0 = AS_NUMERIC(t0));
  PROTECT(s = AS_NUMERIC(s));
  PROTECT(N = AS_INTEGER(N));
  
 // T = *REAL(t0);
  S = *REAL(s);
  n = *INTEGER(N);
  m = *INTEGER(M);
  DELTA = *REAL(delta);

  if(m>1)
   PROTECT(X = allocMatrix(REALSXP, n+1, m));
  else
   PROTECT(X = NEW_NUMERIC(n+1));
  rx0 = REAL(x0);  
  rX = REAL(X);
  for(j=0; j<m; j++)
   rX[j*(n+1)] = rx0[j]; 

  GetRNGstate();
  for(i=1; i< n+1; i++){
   for(j=0; j<m;j++){    
	tmp = rX[i + j*(n+1) - 1];  
    D = feval(DELTA,tmp,d,rho);
	Dx = feval(DELTA,tmp,dx,rho);
    Dxx = feval(DELTA,tmp,dxx,rho);
    Dt = feval(DELTA,tmp,dt,rho);
	Mx = S*S * Dxx/2.0 + Dt;
	Ex = tmp + D*(exp(Dx*DELTA)-1.0)/Dx + 
	       Mx*(exp(Dx*DELTA) -1.0 -Dx*DELTA)/(Dx*Dx); 
    Vx = S*sqrt((exp(2.0*Dx*DELTA)-1.0)/(2.0*Dx));
	rX[i + j*(n+1)] = rnorm(Ex,Vx);
	}
  }   
  PutRNGstate();

  UNPROTECT(6);
  return(X);
}

/* approximate transition densities */

/* Hermite polynomials up to order 6*/
double H0(double z){ 
 return(1.0);
}
double H1(double z){ 
 return(-z);
}
double H2(double z){ 
 return(-1.0+z*z);
}
double H3(double z){ 
 return(3.0*z-z*z*z);
}
double H4(double z){ 
 double z2 = z*z;
 return(3.0-6.0*z2+z2*z2);
} 
double H5(double z){
 double z2=z*z;
 return(-15.0*z+10.0*z*z2-z*z2*z2);
}
double H6(double z){
 double z2=z*z;
 return(-15.0+45.0*z2-15.0*z2*z2+z2*z2*z2);
}




/* HPloglik: likelihood approximation by Hermite poynomials 
   Ait-Sahalia, Y. (2002) "Maximum likelihood estimation of
   discretely sampled diffusions: a clsed form approximation 
   approach", Econometrica, 70, 223-262

	Delta: time step
	X: data
	theta: vector of parameters
	M0: drift of the trans. diffusion
	M1,..., M6: derivatives of M0 up to order 6
	F: transform function
	S: diffusion coefficient
	F, S, M0-M6 are functions of (x,theta)
    There is no sanity check on these functions. Be warned!!!
    see text, Chapter 3.
	
	To obtain the transition density, just call HP with X
	of length 2.
	
	function returns the log-likelihood
*/   
   

SEXP HPloglik(SEXP delta, SEXP X, SEXP theta, SEXP M0, SEXP M1,
                  SEXP M2, SEXP M3, SEXP M4, SEXP M5, SEXP M6, 
				  SEXP F, SEXP S, SEXP rho){
 double f, y0, Delta;
 double sd, ssd;
 double mu0, mu1, mu2, mu3, mu4, mu5, mu6;
 double val=0.0;
 int i, n;
 double *x;
 SEXP ans;
 
 if(!isNumeric(X)) error("`X' must be numeric");
 if(!isNumeric(delta)) error("`delta' must be numeric");

 PROTECT(ans = NEW_NUMERIC(1));
 PROTECT(delta = AS_NUMERIC(delta));
 PROTECT(X = AS_NUMERIC(X));
 PROTECT(theta);
 PROTECT(M0);
 PROTECT(M1);
 PROTECT(M2);
 PROTECT(M3);
 PROTECT(M4);
 PROTECT(M5);
 PROTECT(M6);
 PROTECT(F);
 PROTECT(S);
 PROTECT(theta);
 
 Delta = *REAL(delta);
 x = REAL(X);
 n = length(X);
 sd = sqrt(Delta);

 for(i=1; i<n; i++){   
  y0 = ftheta(0,x[i-1],theta, F, rho);
  f = ftheta(0,x[i],theta, F, rho);
  ssd = ftheta(0,x[i], theta, S, rho)*sd;
  mu0 = ftheta(0,y0, theta, M0, rho);
  mu1 = ftheta(0,y0, theta, M1, rho);  
  mu2 = ftheta(0,y0, theta, M2, rho);  
  mu3 = ftheta(0,y0, theta, M3, rho);  
  mu4 = ftheta(0,y0, theta, M4, rho);  
  mu5 = ftheta(0,y0, theta, M5, rho);  
  mu6 = ftheta(0,y0, theta, M6, rho);  
  val+= log(cHP(Delta, mu0, mu1, mu2, mu3, mu4, mu5, mu6, (f-y0)/sd, ssd));
 }
  
 REAL(ans)[0] = val; 
 UNPROTECT(14);
 return(ans);
}

SEXP EULERloglik(SEXP delta, SEXP X, SEXP theta, SEXP d, SEXP s, SEXP rho){
 double Delta, sd;
 double val=0.0;
 int i, n;
 double *x;
 SEXP ans;
 
 if(!isNumeric(X)) error("`X' must be numeric");
 if(!isNumeric(delta)) error("`delta' must be numeric");

 PROTECT(ans = NEW_NUMERIC(1));
 PROTECT(delta = AS_NUMERIC(delta));
 PROTECT(X = AS_NUMERIC(X));
 PROTECT(theta);
 PROTECT(d);
 PROTECT(s);
 PROTECT(theta);
 
 Delta = *REAL(delta);
 x = REAL(X);
 n = length(X);
 sd = sqrt(Delta);

 for(i=1; i<n; i++){  
  val += dnorm(x[i], x[i-1] + ftheta(0,x[i-1], theta, d, rho)*Delta, 
                sd* ftheta(0, x[i-1], theta, s, rho), TRUE);  
 }
  
 REAL(ans)[0] = val; 
 UNPROTECT(7);
 return(ans);
}


SEXP dcSim(SEXP X, SEXP Y, SEXP delta, SEXP d, SEXP s, SEXP theta, SEXP N, SEXP M, SEXP rho){
 double Delta, sd, tmp, tmp1, x1, x2, z;
 int i, k, m;
 double x, y;
 int NN, MM;
 SEXP ans;
 
 if(!isNumeric(X)) error("`X' must be numeric");
 if(!isNumeric(Y)) error("`Y' must be numeric");
 if(!isNumeric(delta)) error("`delta' must be numeric");
 if(!isInteger(N)) error("`N' must be integer");
 if(!isInteger(M)) error("`M' must be integer");

 PROTECT(delta = AS_NUMERIC(delta));
 PROTECT(X = AS_NUMERIC(X));
 PROTECT(Y = AS_NUMERIC(Y));
 PROTECT(N = AS_INTEGER(N));
 PROTECT(M = AS_INTEGER(M));
 PROTECT(theta);
 PROTECT(d);
 PROTECT(s);
 PROTECT(theta);
 
 NN = *INTEGER(N);
 Delta = *REAL(delta)/(double)NN;
 NN--;
 MM = *INTEGER(M);
 x = *REAL(X);
 y = *REAL(Y);
 PROTECT(ans = NEW_NUMERIC(1));
 sd = sqrt(Delta);

 GetRNGstate();
 tmp = 0.0;
 k = 0;

 for(m=0; m<MM-1; m+=2){ /* MC iterations */
  x1 = x2 = x;
  for(i=1; i<NN; i++){  /* number of sub-intervals */
   z = rnorm(0,1);
   x1 += ftheta(0, x1, theta, d, rho)*Delta + ftheta(0, x1, theta, s, rho)*sd*z; 
   x2 += ftheta(0, x2, theta, d, rho)*Delta - ftheta(0, x2, theta, s, rho)*sd*z;    
  }
  tmp1 = dnorm(y, x1 + ftheta(0, x1, theta, d, rho)*Delta,
                sd*ftheta(0, x1, theta, s, rho), FALSE);
  if(!isnan(tmp1)){				
   tmp += tmp1;
   k++;
  }
  tmp1 = dnorm(y, x2 + ftheta(0, x2, theta, d, rho)*Delta,
                sd*ftheta(0, x2, theta, s, rho), FALSE);
  if(!isnan(tmp1)){				
   tmp += tmp1;
   k++;
   }
 } /* MM */
 REAL(ans)[0] = tmp/k;
 
 PutRNGstate();
  
 UNPROTECT(10);
 return(ans);
}



SEXP SIMloglik(SEXP X, SEXP delta, SEXP d, SEXP s, SEXP theta, SEXP N, SEXP M, SEXP rho){
 double Delta, sd, tmp, tmp1, x1, x2, z;
 int h, i, k, m, n;
 double *x, val;
 int NN, MM;
 SEXP ans;
 
 if(!isNumeric(X)) error("`X' must be numeric");
 if(!isNumeric(delta)) error("`delta' must be numeric");
 if(!isInteger(N)) error("`N' must be integer");
 if(!isInteger(M)) error("`M' must be integer");

 PROTECT(delta = AS_NUMERIC(delta));
 PROTECT(X = AS_NUMERIC(X));
 PROTECT(N = AS_INTEGER(N));
 PROTECT(M = AS_INTEGER(M));
 PROTECT(theta);
 PROTECT(d);
 PROTECT(s);
 PROTECT(theta);
 n = length(X);
  
 NN = *INTEGER(N);
 Delta = *REAL(delta)/(double)NN;
 NN--;
 MM = *INTEGER(M);
 x = REAL(X);
 PROTECT(ans = NEW_NUMERIC(1));
 sd = sqrt(Delta);

 GetRNGstate();
 val = 0;
 for(h=1; h<n; h++){
  tmp = 0.0;
  k = 0;

  for(m=0; m<MM-1; m+=2){ /* MC iterations */
   x1 = x2 = x[h-1];
   for(i=1; i<NN; i++){  /* number of sub-intervals */
    z = rnorm(0,1);
    x1 += ftheta(0, x1, theta, d, rho)*Delta + ftheta(0, x1, theta, s, rho)*sd*z; 
    x2 += ftheta(0, x2, theta, d, rho)*Delta - ftheta(0, x2, theta, s, rho)*sd*z;    
   }
   tmp1 = dnorm(x[h], x1 + ftheta(0, x1, theta, d, rho)*Delta,
                sd*ftheta(0, x1, theta, s, rho), FALSE);
   if(!isnan(tmp1)){				
    tmp += tmp1;
    k++;
   }
   tmp1 = dnorm(x[h], x2 + ftheta(0, x2, theta, d, rho)*Delta,
                sd*ftheta(0, x2, theta, s, rho), FALSE);
   if(!isnan(tmp1)){				
    tmp += tmp1;
    k++;
    }
  } /* MM */
  val += log(tmp/k);
 }
 REAL(ans)[0] = val;
 PutRNGstate();
  
 UNPROTECT(9);
 return(ans);
}


double cHP(double Delta, double mu0, double mu1, double mu2, double mu3,
						   double mu4, double mu5, double mu6, 
						   double z, double ssd){
 double eta1, eta2, eta3, eta4, eta5, eta6;
 double mu02, mu03, mu04, mu05, mu06;
 double mu12, mu13, mu22;
 double val;
 
 mu02 = pow(mu0,2.0);
 mu03 = pow(mu0,3.0);
 mu04 = pow(mu0,4.0);
 mu05 = pow(mu0,5.0);
 mu06 = pow(mu0,6.0);
 mu12 = pow(mu1,2.0);
 mu13 = pow(mu1,3.0);
 mu22 = pow(mu2,2.0);
 
 eta1 = -mu0*sqrt(Delta) - (2.0*mu0*mu1+mu2)*pow(Delta, 1.5)/4.0
  -(4.0*mu0*mu12 + 4.0*mu02*mu2 + 6.0*mu1*mu2 + 
    4.0*mu0*mu3 + mu4)*pow(Delta,2.5)/24.0;

 eta2 = (mu02+mu1)*Delta/2.0 + (6.0*mu02*mu1 + 4.0*mu12
  + 7.0*mu0*mu2 + 2*mu3)*pow(Delta,2.0)/12.0 +(28.0*mu02*mu12
  + 28.0*mu02*mu3 + 16.0*mu13 
  + 16.0*mu03*mu2 + 88.0*mu0*mu1*mu2 + 21.0*mu22 + 32.0*mu1*mu3 
  + 16.0*mu0*mu4 + 3.0*mu5)*pow(Delta,3.0)/96.0;

 eta3 = -(mu03 + 3.0*mu0*mu1 + mu2)*pow(Delta,1.5)/6.0 - (12.0*mu03 * mu1 
  + 28.0*mu0*mu12 + 22.0*mu02*mu2 + 24.0*mu1*mu2 + 14.0*mu0*mu3 
  + 3.0*mu4)*pow(Delta,2.5)/48.0;

 eta4 = (mu04 + 6.0*mu02 * mu1 + 3.0*mu12 + 4.0*mu0*mu2 + mu3)*pow(Delta,2.0)/24.0
  + (20.0*mu04*mu1 + 50.0*mu03*mu2 + 100.0*mu02*mu12 + 50.0*mu02*mu3 
  + 23.0*mu0*mu4 + 180.0*mu0*mu1*mu2 + 40.0*mu13 + 34.0*mu22 + 52.0*mu1*mu3
  + 4.0*mu5)*pow(Delta,3.0)/240.0;

 eta5 = -(mu05 + 10.0*mu03*mu1 + 15.0*mu0*mu12 + 10.0*mu02*mu2
  + 10.0*mu1*mu2 + 5.0*mu0*mu3 + mu4)*pow(Delta,2.5)/120.0;
 
 eta6 = (mu06 + 15.0*mu04*mu1 + 15.0*mu13 + 20.0*mu03*mu2 + 15.0*mu0*mu3 
  + 45.0*mu02*mu12 + 10.0*mu22 + 15.0*mu02*mu3 + 60.0*mu0*mu1*mu2
  + 6.0*mu0*mu4 + mu5)*pow(Delta,3.0)/720.0;

 val = dnorm(z, 0.0, 1.0, FALSE) * (1.0 + eta1*H1(z) + eta2*H2(z) + eta3*H3(z) 
        + eta4*H4(z) + eta5*H5(z) + eta6*H6(z))/ssd;

 return(val);
}



static R_CMethodDef R_CDef[] = {
   {"sde_sim_euler", (DL_FUNC)&sde_sim_euler, 12},
   {"sde_sim_milstein", (DL_FUNC)&sde_sim_milstein, 9},
   {"sde_sim_milstein2", (DL_FUNC)&sde_sim_milstein2, 12},
   {"sde_sim_KPS", (DL_FUNC)&sde_sim_KPS, 14},
   {"sde_sim_cdist", (DL_FUNC)&sde_sim_cdist, 7},
   {"sde_sim_ozaki", (DL_FUNC)&sde_sim_ozaki, 9},
   {"sde_sim_shoji", (DL_FUNC)&sde_sim_shoji, 11},
   {"HPloglik", (DL_FUNC)&HPloglik, 13},
   {"EULERloglik", (DL_FUNC)&EULERloglik, 6},
   {"SIMloglik", (DL_FUNC)&SIMloglik, 8},
   {"dcSim", (DL_FUNC)&dcSim, 9},
   {NULL, NULL, 0},
};

void
R_init_ifs(DllInfo *info)
{
    R_registerRoutines(info, R_CDef, NULL, NULL, NULL);
    R_useDynamicSymbols(info, TRUE);
}


/* Accessor functions. For internal use only. */

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
    val = *REAL(eval(R_fcall, rho));
    UNPROTECT(3);
		
    return(val);

}

/* ftheta: this is intended to calculate   f(t, x,theta) where
   't' and 'x' are C doubles and 'theta' is a SEXP coming from the
   R workspace.
*/   

double ftheta(double t, double x, SEXP theta, SEXP f, SEXP rho)
{
    double val= 0.0;
    SEXP R_fcall, tpar, xpar; 
	
	PROTECT(theta);
	PROTECT(tpar = allocVector(REALSXP, 1));
	PROTECT(xpar = allocVector(REALSXP, 1));
    REAL(tpar)[0] = t;
    REAL(xpar)[0] = x;
   
	PROTECT(R_fcall = allocList(4));
	SETCAR(R_fcall, f);
	SET_TYPEOF(R_fcall, LANGSXP);

	SETCADR(R_fcall, tpar);
	SETCADDR(R_fcall, xpar);
	SETCADDDR(R_fcall, theta);
    val = *REAL(eval(R_fcall, rho));
    UNPROTECT(4);
		
    return(val);
}

