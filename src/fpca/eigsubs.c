#include <stdio.h>
#include <limits.h>
#include <math.h>  
#include <nicklib.h> 
#include "eigsubs.h" 
/* ********************************************************************* */

void packsym(double *pmat, double *mat, int n) ;


void eigvals(double *mat, double *evals, int n) 
{
	 double *pmat ;  

	 ZALLOC(pmat, (n*n), double) ;
	 vst(mat, mat, -1.0, n*n) ;
	 packsym(pmat, mat, n) ;
	 

         eigx_(pmat, evals, &n) ;
	 free(pmat) ;
	 vst(mat, mat, -1.0, n*n) ;
	 vst(evals, evals, -1.0, n) ;
}
void eigvecs(double *mat, double *evals, double *evecs, int n) 
{
	 double *pmat ;  

	 ZALLOC(pmat, (n*n), double) ;
	 vst(mat, mat, -1.0, n*n) ;
	 packsym(pmat, mat, n) ;

         eigxv_(pmat, evals, evecs, &n) ;
	 free(pmat) ;
	 vst(mat, mat, -1.0, n*n) ;
	 vst(evals, evals, -1.0, n) ;
}
void eigb(double *lam, double *a, double *b, int n) 
// bidiagonal matrix  
{  
   double *w, *d, *e, *ww ;

//1: form B B'  
   ZALLOC(w, 3*n, double) ; 
   d = w ; 
   e = w+n ;  
   ww = e+n ;
   vvt(e, a, b, n-1) ;  
   vvt(d, a, a, n) ;
   vvt(ww, b, b, n-1) ;
   vvp(d+1, d+1, ww, n) ;
//2: call tridiag solver
   eigc(lam, d, e, n) ;

   free(w) ;
}

void eigc(double *lam, double *a, double *b, int n) 
// tridiagonal matrix  
{  

   double *w, *d, *e ;  
   int nn = n, info ;
   ZALLOC(w, 2*n, double) ; 
   d = w ; 
   e = w+n ;  
   copyarr(a, d, n) ;
   copyarr(b, e, n-1) ;
   vst(w, w, -1.0, 2*n) ;
   dsterf_(&nn, d, e, &info) ;
   if (info != 0) fatalx("(eigc) %d\n", info) ;
   vst(lam, d, -1.0, n) ;

   free(w) ;
}

   

   
void
packsym(double *pmat, double *mat, int n) 
	//  lapack L mode (fortran)
{ 
	int i, j, k = 0 ;
	for (i=0; i<n; i++)  {  
          for (j=i; j<n; j++) { 
             pmat[k] = mat[i*n+j] ;
	     ++k ;
	  }
	}
}

