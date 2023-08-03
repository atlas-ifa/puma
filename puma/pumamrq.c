#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXPAR  32		/* Max number of parameters */
#define MAXITER 20		/* Max iterations */
#define QFRAC   1e-8		/* Fractional change for quitting */
// #define QFRAC   1e-12		/* Fractional change for quitting */

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

/* REAL is 16 byte floating point */
typedef long double REAL;
#define SQRT(x) (sqrtl(x))
#define POW(x,y) (powl(x,y))
#define FABS(x) (fabsl(x))

/* Uncomment only one for numerical derivatives */
//#define OLDDERIV	/* Compute cross derivatives in ++ quadrant */
//#define SYMDERIV	/* Compute cross derivatives symmetrically */
#define PDERIV		/* Compute cross derivatives from parabolic fit */

/* Invert a matrix */
int invert(int n, REAL *a, REAL *det);


/* Definition of function that mini() minimizes */
typedef int MINIFUNC(int npar, REAL *par, int *usepar, REAL *value,
			 int doderiv, REAL *deriv, REAL *curve);

/* Black box minimization using Levenberg-Marquard algorithm */
int mini(int npar, REAL *par, int *usepar,
	 MINIFUNC func, REAL *cov, int nderiv, int *lambda, int verbose);

/* Debug print of status and variables */
int vprint(int n, REAL *a, REAL f, int niter, int lambda);

/* Show the covariance matrix */
int covprint(int npar, REAL *cov);

/* Black box minimization using Levenberg-Marquard algorithm */
int mini(int npar, REAL *par, int *usepar,
	 MINIFUNC func, REAL *cov, int nderiv, int *lambda, int verbose)
{
/*
*     This is John Tonry's black box minimization (and fitting) program,
*     implementing a Marquard algorithm.
*     Transcribed to C, 020904.
*     Revision 2.0, 11/17/82.
*
*     MINI's arguments are as follows:
*     NPAR - The number of parameters to be varied in searching for a
*        minimum.
*     PAR - A double vector, it serves three purposes:
*        (1) It passes MINI the initial estimate for the parameters
*        (2) It passes arguments to the function to be minimized
*        (3) It returns the minimum parameter values
*     USEPAR - set to 0/1 to turn on that parameter for variation
*     FUNC - The function to be minimized. FUNC takes NPAR, PAR and 
*        USEPAR as arguments and return the value of the function as
*        VALUE.  If DODERIV = 0, this is all that is wanted.  If
*        DODERIV = 1, the first partial derivatives wrt the active
*        parameters should be returned in DERIV (close packed).  If
*        DODERIV = 2, the second partial derivatives wrt the active
*        parameters should be returned in CURVE (close packed).
*     COV - A NxN real*8 matrix in which the covariance matrix of the
*        fit is returned.
*     NDERIV - Variable to govern how MINI gets its derivatives:
*        NDERIV = -N numerical derivatives, derivative step 10^-N
*        NDERIV = 0 for a function of arbitrary A and numerical derivatives.
*        NDERIV = 1 for a function of arbitrary A, analytic first derivatives
*                  provided by DFUNK, and numerical second derivatives.
*        NDERIV = 2 for a function of arbitrary A, analytic first and second
*                  derivatives provided by DFUNK and by D2FUNK.
*        NDERIV = 3 for a function that is quadratic in A, and MINI
*                  will iterate once, computing the minimum exactly.
*     VERBOSE - governs whether MINI will print each iteration
*
*     LAMBDA - starting value for lambda (suggest -2), returns final value
*
*     Return value - 0 successful completion
*                   -1 error from invert: singular matrix
*                   -2 error from invert: matrix too large
*                   -3 too many parameters
*                    N maximum iteration reached (N)
*
*     Descriptions of some of the variables:
*     PAR - Argument for the function
*     A0 - Current guess for the minimum
*     AI - increments for a0 in computing derivatives
*     DA - Vector from A0 to new guess for the minimum
*     DF - First derivatives of the function
*     D2F - Second derivatives of the function
*     LAMBDA - Governs mix of gradient and analytic searches
*     ITER - Maximum number of iterations
*     QFRAC - Maximum fractional change for successful exit
*
*     The calling program should be as follows (eg):
*********************************************************
*     REAL*8 A(4)
*     REAL*8 COV(4,4)
*     EXTERNAL CHI, DCHI, D2CHI
*     ... (Initialize A to the guess for the minimum)
*     CALL MINI(4,A,CHI,DCHI,D2CHI,COV,NDERIV)
*     ...
*     FUNCTION CHI(A)
*     REAL*8 A(4), CHI
*     ... (define the function)
*     SUBROUTINE DCHI(A,DF)
*     REAL*8 A(4), DF(4)
*     ... (dummy or compute the derivatives)
*     SUBROUTINE D2CHI(A,COV)
*     REAL*8 A(4), COV(4,4)
*     ... (dummy or compute the derivatives)
*************************************************************
*/
   REAL df[MAXPAR], a0[MAXPAR], da[MAXPAR], ai[MAXPAR];
   REAL d2f[MAXPAR*MAXPAR], fnow, fthen, fminus, curve, det, err, base=10.0;
#ifdef SYMDERIV
   REAL fpm, fmp, fmm;
#endif
// OLD versions
//   double dfrac=0.02, places=1e-7, vary=1e-5;
//   double dfrac=0.02, places=1e-6, vary=1e-4;
// Current working...  Pretty efficient
   REAL dfrac=0.01, places=1e-7, vary=1e-5;
   int on[MAXPAR];
   int i, j, k, lam, nuse, error, iter=0, lamit;
      
   if(npar > MAXPAR) {
      fprintf(stderr, "Too many parameters, max = %d\n", MAXPAR);
      return(-3);
   }

   if(nderiv < 0) {
      dfrac = places = vary = pow(10.0, (double)nderiv);
      nderiv = 0;
   }

/*     Define a few parameters */
   lam = *lambda;
   for(i=0, nuse=0; i<npar; i++) {
      if(usepar[i]) {
	 on[nuse] = i;
	 nuse++;
      }
   }

/*     If NDERIV = 3, compute the minimum directly and exit. */
   if(nderiv == 3) {
      for(i=0; i<nuse; i++) par[on[i]] = 0.0;
      error = func(npar, par, usepar, &fnow, 0, NULL, NULL);
      for(i=0; i<nuse; i++) {
	 par[on[i]] = 1.0;
	 error = func(npar, par, usepar, &df[i], 0, NULL, NULL);
	 for(j=0; j<=i; j++) {
	    par[on[j]] = par[on[j]] + 1;
	    error = func(npar, par, usepar, &d2f[i+nuse*j], 0, NULL, NULL);
	    d2f[i+nuse*j] += fnow - df[i] - df[j];
	    cov[i+nuse*j] = cov[j+nuse*i] = d2f[j+nuse*i] = d2f[i+nuse*j];
	    par[on[j]] = 0.0;
	 }
	 par[on[i]] = 0.0;
      }
      for(i=0; i<nuse; i++) df[i] -= 0.5*d2f[i+i*nuse] + fnow;

      if( (error=invert(nuse, cov, &det)) != 0) return(error);

      for(i=0; i<nuse; i++) {
	 for(j=0, par[on[i]]=0.0; j<nuse; j++) 
	    par[on[i]] -= cov[i+nuse*j] * df[j];
      }
      error = func(npar, par, usepar, &fnow, 0, NULL, NULL);

      goto DOCOVAR;
      fprintf(stderr, "minimrq: Should never see this...\n");
   }

/*     Initialize A0 */
   for(i=0; i<npar; i++) a0[i] = par[i];
   error = func(npar, par, usepar, &fnow, 0, NULL, NULL);

   if(verbose) vprint(npar, a0, fnow, 0, lam);

/*     Initialize AI */
   for(i=0; i<nuse; i++) {
      ai[on[i]] = ABS(vary*a0[on[i]]);
      if(ai[on[i]] == 0.0) ai[on[i]] = vary;
   }

/*     Begin iteration to find minimum */
   for(iter=0; iter<MAXITER; iter++) {

      fthen = fnow;

/* Compute the function derivatives. */
      for(j=0; j<nuse; j++) par[on[j]] = a0[on[j]];

      switch(nderiv) {

/* First case: NDERIV = 0 so entirely numerical derivatives are required */
	 case 0:
/* First the 1st derivatives */
	    for(j=0; j<nuse; j++) {
	       par[on[j]] = a0[on[j]] + ai[on[j]];
	       error = func(npar, par, usepar, &df[j], 0, NULL, NULL);
	       par[on[j]] = a0[on[j]];
	    }

/* The off-diagonal 2nd derivatives */
#ifndef PDERIV
	    for(j=1; j<nuse; j++) {
	       for(k=0; k<j; k++) {
		  par[on[k]] = a0[on[k]] + ai[on[k]];
		  par[on[j]] = a0[on[j]] + ai[on[j]];
		  error = func(npar, par, usepar, &d2f[j+k*nuse], 0, NULL, NULL);
#ifdef OLDDERIV
// Use just the ++ quadrant instead of symmetric around the origin
		  d2f[j+k*nuse] = d2f[k+j*nuse] = 
		     (d2f[j+k*nuse]-df[k]-df[j]+fnow) / (ai[on[j]]*ai[on[k]]);
#endif
#ifdef SYMDERIV
// Use the full area around the origin
		  par[on[j]] = a0[on[j]] - ai[on[j]];
		  error = func(npar, par, usepar, &fpm, 0, NULL, NULL);
		  par[on[k]] = a0[on[k]] - ai[on[k]];
		  error = func(npar, par, usepar, &fmm, 0, NULL, NULL);
		  par[on[k]] = a0[on[k]] - ai[on[k]];
		  par[on[j]] = a0[on[j]] + ai[on[j]];
		  error = func(npar, par, usepar, &fmp, 0, NULL, NULL);
		  d2f[j+k*nuse] = d2f[k+j*nuse] = 
		     (d2f[j+k*nuse]+fmm-fmp-fpm) / (4*ai[on[j]]*ai[on[k]]);
#endif		  
		  par[on[k]] = a0[on[k]];
		  par[on[j]] = a0[on[j]];
	       }
	    }
#endif

/* The on-diagonal 2nd derivatives, and fix the 1st ones. */
	    for(j=0; j<nuse; j++) {
	       par[on[j]] = a0[on[j]] - ai[on[j]];
	       error = func(npar, par, usepar, &fminus, 0, NULL, NULL);
	       d2f[j+j*nuse] = (fminus+df[j]-2*fnow) / (ai[on[j]]*ai[on[j]]);
	       df[j] = (df[j] - fminus) / (2*ai[on[j]]);
	       par[on[j]] = a0[on[j]];
	    }

/* The off-diagonal 2nd derivatives */
#ifdef PDERIV
	    for(j=1; j<nuse; j++) {
	       for(k=0; k<j; k++) {
		  par[on[k]] = a0[on[k]] + ai[on[k]];
		  par[on[j]] = a0[on[j]] + ai[on[j]];
		  error = func(npar, par, usepar, &d2f[j+k*nuse], 0, NULL, NULL);
		  d2f[j+k*nuse] = d2f[k+j*nuse] = 
		     ((d2f[j+k*nuse] - fnow)/(ai[on[j]]*ai[on[k]]) - 
		      df[j]/ai[on[k]] - df[k]/ai[on[j]] -
		      0.5*d2f[j+j*nuse]*ai[on[j]]/ai[on[k]] -
		      0.5*d2f[k+k*nuse]*ai[on[k]]/ai[on[j]]) ;
		  par[on[k]] = a0[on[k]];
		  par[on[j]] = a0[on[j]];
	       }
	    }
#endif
	    break;

/* Second case: NDERIV = 1 so analytic first derivatives are available */
	 case 1:
	    error = func(npar, par, usepar, &fminus, 1, df, NULL);
	    for(j=0; j<nuse; j++) {
	       par[on[j]] = a0[on[j]] + ai[on[j]];
	       error = func(npar, par, usepar, &fminus, 1, da, NULL);
	       par[on[j]] = a0[on[j]];
	       for(i=0; i<=j; i++) {
		  d2f[i+j*nuse] = d2f[j+i*nuse] = (da[i]-df[i]) / ai[on[j]];
	       }
	    }
	    break;

/* Third case: NDERIV = 2 so analytic derivatives are available */
	 case 2:
	    error = func(npar, par, usepar, &fminus, 2, df, d2f);
	    break;

	 default:
	    fprintf(stderr, "unknown nderiv = %d\n", nderiv);
	    return(-3);
      }

/* Compute better estimates for the increments. */
      for(j=0; j<nuse; j++) {
	 curve = d2f[j+j*nuse];
	 if(curve == 0.0) curve = vary;
	 ai[on[j]] = SQRT(POW(df[j]*dfrac/curve,2.0)+ABS(fnow*places/curve));
//	 printf(" %12.4Le", ai[on[j]]);
      }
//      printf("\n");
/*     Begin loop to find a direction along which function decreases */
      for(lamit=0; lamit<MAXITER; lamit++) {

/*     Get weight matrix */
	 for(j=0; j<nuse; j++) {
	    for(i=0; i<j; i++) cov[i+j*nuse] = cov[j+i*nuse] = d2f[i+j*nuse];
	    cov[j+j*nuse] = ABS(d2f[j+j*nuse]*(1+POW(base,(REAL)lam)));
	 }
	 if( (error=invert(nuse, cov, &det)) != 0) return(error);

/*     Multiply to get dA */
	 for(j=0; j<nuse; j++) {
	    for(i=0,da[j]=0.0; i<nuse; i++) da[j] -= cov[j+i*nuse]*df[i];
	 }
/*     Now get new function value */
	 for(j=0; j<nuse; j++) par[on[j]] = a0[on[j]] + da[j];
	 error = func(npar, par, usepar, &fnow, 0, NULL, NULL);

//	 for(j=0; j<nuse; j++) printf(" %12.4e", df[j]);
//	 printf("\n");
//	 vprint(npar, par, fnow, lamit, lam);
/*
 *     Test for whether the function has decreased
 *     If so, adopt the new point and decrement LAMBDA
 *     Else, increment LAMBDA, and get a new weight matrix
 */
	 if(fnow < fthen) break;
	 lam++;
      }
/*     Normal exit, the function at A0 + DA is less than at A0 */
      for(j=0; j<nuse; j++) a0[on[j]] = par[on[j]];
      lam--;
/*
 *     Print the current status and test to see whether the function
 *     has varied fractionally less than QFRAC.
 */
      if(verbose) vprint(npar, a0, fnow, iter, lam);
      if(ABS(fthen-fnow)/fnow < QFRAC) break;
   }
/*     This is the final computation of the covariance matrix */
/*     Quit if no minimum was found in the allowed number of iterations */

DOCOVAR:

/*     Finally, compute the covariance matrix */
   for(j=0; j<nuse; j++) {
      for(k=0; k<nuse; k++) cov[k+j*nuse] = d2f[k+j*nuse] / 2;
   }

   if(verbose > 1) {
      printf("D2F:\n");
      covprint(nuse, cov);
   }

   if( (error=invert(nuse, cov, &det)) != 0) return(error);
   for(j=0; j<nuse; j++) {
      err = SQRT(ABS(cov[j+j*nuse]));
      if(cov[j+j*nuse] < 0) err = -err;
      cov[j+j*nuse] = err;
   }

   for(j=1; j<nuse; j++) {
      for(k=0; k<j; k++) {
	 cov[k+j*nuse] /= cov[k+k*nuse]*cov[j+j*nuse];
	 cov[j+k*nuse] = cov[k+j*nuse];
      }
   }
   if(verbose > 1) {
      printf("cov:\n");
      covprint(nuse, cov);
   }

   *lambda = lam;

   if(iter >= MAXITER) {
      if(verbose) printf("Maximum iteration exceeded\n");
      return(iter);
   }
   return(0);
}

int vprint(int n, REAL *a, REAL f, int niter, int lambda)
{
   int i;
   for(i=0; i<n; i++) {
      if((i%7) == 0) printf(" A(I) =");
      printf("%11.4Lg", a[i]);
      if((i%7) == 6) printf("\n");
   }
   if(((i-1)%7) != 6) printf("\n");
   printf(" F = %14.7Lg     ITER = %3d    LAMBDA = %3d\n", f, niter, lambda);
   return(0);
}

int covprint(int npar, REAL *cov)
{
   int i, j;
   for(j=0; j<npar; j++) {
      for(i=0; i<npar; i++) 
	 printf("%9.4Lf%c", cov[i+j*npar], (i==npar-1) ? '\n' : ' ');
   }
   return(0);
}


#define MAXSIZE 1000
#define MAXDEC 1e30	/* Maximum decrement in pivot value for singular */

int invert(int n, REAL *a, REAL *det)
{

/*
      SUBROUTINE INVERT(N,A,RST,DET)

     Subroutine to invert a matrix, and compute the determinant.
     John Tonry, 9/2/80; transcribed to C 020904.

     INVERT's arguments:

     N - The dimension of the matrix
     A - The NxN matrix to be inverted. Upon successful inversion, A
         contains the inverse. A must be real*8
     DET - The determinant of the matrix. DET is set to 0 for a
         singular matrix, and in that case, A contains garbage.

     return values are 0 -- normal return
                      -2 -- matrix too large
                      -1 -- singular matrix
*/
   REAL save, pivot, onrow, cprev, cnow, decr;
   short int rst[MAXSIZE][2];
   int i, j, k, l, mrank, isign, nrow=0, ncol=0;

   if(n > MAXSIZE) {
      fprintf(stderr, "Matrix size too big: %d > %d\n", n, MAXSIZE);
      return(-2);
   }
   mrank = 0;
   isign = 1;
   *det = 0.0;
   for(j=0; j<n; j++) rst[j][0] = rst[j][1] = -1;
/* Loop over columns, reducing each */
   for(i=0; i<n; i++) {

/* Find the pivot element */
      pivot = -1.0;
      for(j=0; j<n; j++) {
	 if(rst[j][0] != -1) continue;
	 for(k=0; k<n; k++) {
	    if(rst[k][0] != -1) continue;
	    if(pivot >= ABS(a[j+k*n])) continue;
	    pivot = ABS(a[j+k*n]);
	    nrow = j;
	    ncol = k;
	 }
      }
      pivot = a[nrow+ncol*n];
      if(pivot == 0.0) {
	 *det = 0;
	 return(-1);
      }
      rst[ncol][0] = nrow;
      rst[ncol][1] = i;
/* Swap pivot element onto the diagonal */
      for(k=0; k<n; k++) {
	 save = a[nrow+k*n];
	 a[nrow+k*n] = a[ncol+k*n];
	 a[ncol+k*n] = save;
      }
/*   Reduce pivot column */
      for(j=0; j<n; j++) a[j+ncol*n] = -a[j+ncol*n]/pivot;
      a[ncol+ncol*n] = 1/pivot;

/*   Reduce other columns */
      for(k=0; k<n; k++) {
	 if(k == ncol) continue;
/*     Find maximum of column to check for singularity */
	 cprev = 0;
	 for(j=0; j<n; j++) cprev = MAX(cprev,ABS(a[j+k*n]));
/*     Reduce the column */
	 onrow = a[ncol+k*n];
	 a[ncol+k*n] = 0;
	 for(j=0; j<n; j++) a[j+k*n] = a[j+k*n] + onrow*a[j+ncol*n];
/*     Find the new maximum of the column */
	 cnow = 0;
	 for(j=0; j<n; j++) cnow = MAX(cnow,ABS(a[j+k*n]));

/*     Quit if too many figures accuracy were lost (singular) */
	 decr = cprev / cnow;
	 if(cnow == 0.0 || decr > MAXDEC) {
	    *det = 0;
	    return(-1);
	 }
      }
      *det = *det + log(ABS(pivot));
      if(pivot < 0) isign *= -1;
      mrank++;
   }

/*     Now untangle the mess */
   for(j=0; j<n; j++) {
      for(k=0; k<n; k++) {
	 if(rst[k][1] != (n-1-j)) continue;
	 ncol = rst[k][0];
	 if(ncol == k) break;
	 for(l=0; l<n; l++) {
	    save = a[l+ncol*n];
	    a[l+ncol*n] = a[l+k*n];
	    a[l+k*n] = save;
	 }
	 break;
      }
   }
	 
   if(ABS(*det) < 88) *det = isign * exp(*det);
   return(0);
}









/* Reduce a real, symmetric matrix to tridiagonal form */
void tred2(int N, REAL *A, REAL *D, REAL *E);

/* Tridiagonal QL Implicit: eigenvectors and values of a tridiagonal matrix */
void tqli(int N, REAL *Q, REAL *D, REAL *e);

#define SIGN(a,b) ((b) >= 0.0 ? FABS(a) : -FABS(a))
#define SQR(a) ((a) == 0.0 ? 0.0 : ((a)*(a)))

REAL pythag(REAL a, REAL b)
{
   REAL absa, absb;
   absa = FABS(a);
   absb = FABS(b);
   if (absa > absb) return absa*SQRT(1.0+SQR(absb/absa));
   else return (absb == 0.0 ? 0.0 : absb*SQRT(1.0+SQR(absa/absb)));
}

/* Householder reduction of a real, symmetric N x N matrix A.  On
 * output A is replaced by the orthogonal matrix A, effecting the
 * transformation.  D returns the diagonal elements of the tridiagonal
 * matrix, and E the off-diagonal elements, with E[0]=0.  Several
 * statements can be omitted if only eigenvalues are to be found, in
 * which case A contains no useful information on output.
 */

/* Reduce a real, symmetric matrix to tridiagonal form */
void tred2(int N, REAL *A, REAL *D, REAL *E)
{
   int l, k, j, i;
   REAL scale, hh, h, g, f;

   for(i=N-1; i>=1; i--) {
      l = i-1;
      h = scale = 0.0;
      if(l > 0) {
	 for(k=0; k<=l; k++) scale += FABS(A[k+i*N]);
	 if(scale == 0.0) {	/* Skip transformation */
	    E[i] = A[l+i*N];
	 } else {
	    for(k=0; k<=l; k++) {
	       A[k+i*N] /= scale;	/* Use scaled A's for transformation */
	       h += A[k+i*N]*A[k+i*N];	/* Form sigma in h */
	    }
	    f = A[l+i*N];
	    g = (f >= 0.0 ? -SQRT(h) : SQRT(h));
	    E[i] = scale*g;
	    h -= f*g;		/* h is now NR equation 11.2.4 */
	    A[l+i*N] = f - g;	/* Store u in the ith row of A */
	    f = 0.0;
	    for(j=0; j<=l; j++) {
	       A[i+j*N] = A[j+i*N] / h;	/* Store u/H in ith column of A */
	       g = 0.0;			/* Form an element of A*u in G */
	       for(k=0; k<=j; k++) g += A[k+j*N]*A[k+i*N];
	       for(k=j+1; k<=l; k++) g += A[j+k*N]*A[k+i*N];
	       E[j] = g/h;		/* Form element of p in temp unused */
	       f += E[j]*A[j+i*N];	/* element of E */
	    }
	    hh = f / (h+h);		/* Form K, equation 11.2.11 */
	    for(j=0; j<=l; j++) {	/* Form q and store in E overwriting p */
	       f = A[j+i*N];		/* Note that E(l)=E(i-1) survives */
	       E[j] = g = E[j] - hh*f;
	       for(k=0; k<=j; k++) A[k+j*N] -= (f*E[k]+g*A[k+i*N]); /*Reduce A*/
	    }
	 }
      } else
	 E[i] = A[l+i*N];
      D[i] = h;
   }
   D[0] = 0.0;
   E[0] = 0.0;
   for(i=0; i<N; i++) {	/* Begin accumulation of transformation matrices */
      l = i-1;
      if(D[i]) {	/* Block skipped when i=0 */
	 for(j=0; j<=l; j++) {
	    g=0.0;
	    for(k=0; k<=l; k++) g += A[k+i*N]*A[j+k*N];  /*Use u and u/H saved*/
	    for(k=0; k<=l; k++) A[j+k*N] -= g*A[i+k*N];  /*in A to form p*q */
	 }
      }
      D[i] = A[i+i*N];
      A[i+i*N] = 1.0;	/* Reset row and col of A to identity for next iter */
      for(j=0; j<=l; j++) A[i+j*N] = A[j+i*N]=0.0;
   }
}

/* QL algorithm with implicit shifts, to determine the eigenvalues and
 * eigenvectors of a real, symmetric, tridiagonal matrix, or of a
 * real, symmetric matrix previously reduced by tred2.  D is a vector
 * of length N.  On input its first N elements are the diagonal
 * elements of the tridiagonal matrix.  On output it returns the
 * eigenvalues.  The vector e inputs the subdiagonal elements of
 * the tridiagonal matrix, with e[0] arbitrary.  On output e is
 * destroyed.  When finding only the eigenvalues, several lines may
 * be omitted, as noted.  If the eigenvectors of a tridiagonal matrix
 * are desired the N x N matrix Q is input as the identity matrix.
 *
 * If the eigenvectors of a matrix that has been reduced by tred2 are
 * required, then Q is input as the matrix output by tred2.  In either
 * case the kth column of Q returns the normalized eigenvector
 * corresponding to D[k]:  diagonal_form = Q^T * A * Q
 */

/* Tridiagonal QL Implicit: eigenvectors and values of a tridiagonal matrix */
void tqli(int N, REAL *Q, REAL *D, REAL *e)
{
   int m, l, iter, i, k;
   REAL s, r, p, g, f, dd, c, b;

   for(i=1; i<N; i++) e[i-1] = e[i];
   e[N-1] = 0.0;
   for(l=0; l<N; l++) {
      iter = 0;
      do {
/* Look for a small, sub-diagonal element */
	 for(m=l; m<N-1; m++) {
	    dd = FABS(D[m]) + FABS(D[m+1]);
	    if((REAL)(FABS(e[m])+dd) == dd) break;
	 }
	 if(m != l) {
	    if(iter++ == 30) {
	       fprintf(stderr, "Too many iterations in tqli");
	       exit(1);
	    }
	    g = (D[l+1]-D[l]) / (2.0*e[l]);
	    r = pythag(g, 1.0);
	    g = D[m] - D[l] + e[l]/(g+SIGN(r,g));
	    s = c = 1.0;
	    p = 0.0;
	    for(i=m-1; i>=l; i--) {
	       f = s*e[i];
	       b = c*e[i];
	       e[i+1] = (r=pythag(f,g));
	       if(r == 0.0) {
		  D[i+1] -= p;
		  e[m] = 0.0;
		  break;
	       }
	       s = f/r;
	       c = g/r;
	       g = D[i+1]-p;
	       r = (D[i]-g)*s + 2.0*c*b;
	       D[i+1] = g + (p=s*r);
	       g = c*r-b;
/* Form eigenvectors */
	       for(k=0; k<N; k++) {
		  f = Q[i+1+k*N];
		  Q[i+1+k*N] = s*Q[i+k*N] + c*f;
		  Q[i+k*N] = c*Q[i+k*N] - s*f;
	       }
	    }
	    if(r == 0.0 && i >= l) continue;
	    D[l] -= p;
	    e[l] = g;
	    e[m] = 0.0;
	 }
      } while (m != l);
   }

/* Order eigenvalues and eigenvectors from smallest to biggest */
  for(i=0; i<N-1; i++) {
/* Find the smallest remaining */
    k = i;
    p = D[i];
    for(m=i+1; m<N; m++) {
      if( D[m] < p ) {
        k = m;
        p = D[m];
      }
    }

/* Swap eigenvalues and eigenvectors */
    if( k != i ) {
      D[k] = D[i];
      D[i] = p;
      for(m=0; m<N; m++) {
        g        = Q[i+m*N];
        Q[i+m*N] = Q[k+m*N];
        Q[k+m*N] = g;
      }
    }
  }
}

/* Multiply two NxN matrices: A*B = C */
void Mmult(int N, REAL *A, REAL *B, REAL *C)
{
   int i, j, k;
   for(k=0; k<N; k++) {
      for(j=0; j<N; j++) {
	 C[j+k*N] = 0.0;
	 for(i=0; i<N; i++) C[j+k*N] += A[i+k*N] * B[j+i*N];
      }
   }
}

/* Transpose an NxN matrix B = A^T */
void Mtrans(int N, REAL *A, REAL *B)
{
   int i, j;
   REAL tmp;
   for(j=1; j<N; j++) {
      for(i=0; i<j; i++) {
	 tmp = A[i+j*N];
	 B[i+j*N] = A[j+i*N];
	 B[j+i*N] = tmp;
      }
   }
}

/* Copy an NxN matrix B = A */
void Mcpy(int N, REAL *A, REAL *B)
{
   int i, j;
   for(j=0; j<N; j++) {
      for(i=0; i<N; i++) B[i+j*N] = A[i+j*N];
   }
}

/* Multiply a vector by an NxN matrix V = A * U */
void Mvec(int N, REAL *A, REAL *U, REAL *V)
{
   int i, j;
   for(j=0; j<N; j++) {
      for(i=0, V[j]=0.0; i<N; i++) V[j] += A[j+i*N] * U[i];
   }
}

// #define EIGBLAB		/* Lots of blabby output from d2feig? */

/* Get d2f from reduced cov matrix using eigenvectors, update par and cov */
void d2feig(int npar, REAL *par, int *usepar, MINIFUNC func, 
	    REAL *cov, REAL *d2f, REAL *eigvec, REAL *eigval)
{
   int i, j, k, nuse, on[MAXPAR];
   REAL chim, chip, chi0, DFVAR;
   REAL *dpar, *df, *eigdpar;
      
/* Which parameters are actively varied? */
   for(i=0, nuse=0; i<npar; i++) {
      if(usepar[i]) { on[nuse] = i; nuse++; }
   }

/* Parameter increments (in eigenvector basis) */
   dpar = (REAL *)calloc(nuse, sizeof(REAL));
/* Function gradient (in eigenvector basis) */
   df = (REAL *)calloc(nuse, sizeof(REAL));
/* Parameter offset to minimum (in eigenvector basis) */
   eigdpar = (REAL *)calloc(nuse, sizeof(REAL));

/* Undo the "reduced covariance", invert to restore C^-1 = d2f/2 */
   for(j=1; j<nuse; j++) {
      for(i=0; i<j; i++) eigvec[i+nuse*j] = eigvec[j+nuse*i] = 
			    cov[i+nuse*j] * cov[i+nuse*i]*cov[j+nuse*j];
   }
   for(j=0; j<nuse; j++) {
      eigvec[j+nuse*j] = cov[j+nuse*j] > 0 ? cov[j+nuse*j]*cov[j+nuse*j] :
	                                 -cov[j+nuse*j]*cov[j+nuse*j];
   }
   invert(nuse, eigvec, &chip);

/* Compute the eigenvalues and eigenvectors */
   tred2(nuse, eigvec, eigval, cov);
   tqli(nuse, eigvec, eigval, cov);

#ifdef EIGBLAB
   printf("\n\nd2feig(): \n");
   printf("Eigenvalues: \n");
   for(i=0; i<nuse; i++) printf(" %16.8Le", eigval[i]);
   printf("\n");
   printf("Eigenvectors: \n");
   for(j=0; j<nuse; j++) {
      for(i=0; i<nuse; i++) printf(" %16.8Lf", eigvec[i+j*nuse]);
      printf("\n");
   }
#endif

/* Function value at origin */
   func(npar, par, usepar, &chi0, 0, NULL, NULL);

   if(ABS(chi0) > 1e-6) {
/* Amount of function increase where numerical derivatives are evaluated */
//   DFVAR = 1e-3 * chi0;
      DFVAR = 1e-6 * chi0;
   } else {
//   DFVAR = 0.1;
      DFVAR = 1e-4;
   }

#ifdef EIGBLAB
   printf("%13.6Lf %13.6Lf %13.6Lf  %13.9Lf %13.9Lf %13.9Lf  %13.9Lf\n", 
	  par[on[0]], par[on[1]], par[on[2]], 
	  par[on[3]], par[on[4]], par[on[5]], chi0);
   printf("\n");
#endif

/* Parameter increments (eigenvector basis) */
   for(j=0; j<nuse; j++) dpar[j] = SQRT(DFVAR / ABS(eigval[j]));

#ifdef EIGBLAB
   printf("Increments: \n");
   for(i=0; i<nuse; i++) printf(" %16.8Le", dpar[i]);
   printf("\n");
#endif

/* Evaluate function gradient and second derivatives */
   for(j=0; j<nuse; j++) {
      for(i=0; i<nuse; i++) par[on[i]] += dpar[j] * eigvec[j+i*nuse];
      func(npar, par, usepar, &chip, 0, NULL, NULL);
#ifdef EIGBLAB
   printf("Func plus %d:   %16.8Le", j, chip-chi0);
#endif
      for(i=0; i<nuse; i++) par[on[i]] -= 2 * dpar[j] * eigvec[j+i*nuse];
      func(npar, par, usepar, &chim, 0, NULL, NULL);
#ifdef EIGBLAB
   printf("    minus %d:  %16.8Le\n", j, chim-chi0);
#endif
      for(i=0; i<nuse; i++) par[on[i]] += dpar[j] * eigvec[j+i*nuse];
      d2f[j+j*nuse] = (chip + chim - 2*chi0) / (dpar[j]*dpar[j]);
      df[j] = (chip - chim) / (2*dpar[j]);
   }

/* Evaluate mixed second derivatives */
   for(k=1; k<nuse; k++) {
      for(j=0; j<k; j++) {
	 for(i=0; i<nuse; i++) par[on[i]] += dpar[j]*eigvec[j+i*nuse] +
				             dpar[k]*eigvec[k+i*nuse];
	 func(npar, par, usepar, &chip, 0, NULL, NULL);
#ifdef EIGBLAB
	 printf("Func mix %d %d:  %16.8Le\n", j, k, chip-chi0);
#endif
	 for(i=0; i<nuse; i++) par[on[i]] -= dpar[j]*eigvec[j+i*nuse] +
				             dpar[k]*eigvec[k+i*nuse];
	 d2f[j+k*nuse] = d2f[k+j*nuse] = (chip - chi0 - 
             df[j]*dpar[j] - df[k]*dpar[k] - 0.5*d2f[j+j*nuse]*dpar[j]*dpar[j] -
             0.5*d2f[k+k*nuse]*dpar[k]*dpar[k]) / (dpar[j]*dpar[k]);
      }
   }

#ifdef EIGBLAB
   printf("df in eigenvector basis: \n");
   for(i=0; i<nuse; i++) printf(" %16.8Le", df[i]);
   printf("\n");

   printf("d2f in eigenvector basis: \n");
   for(j=0; j<nuse; j++) {
      for(i=0; i<=j; i++) printf(" %16.8Le", d2f[i+j*nuse]);
      printf("\n");
   }
#endif

/* Tune up the parameters */
// I do seem to find a slightly better func, but somehow iteration gets worse
   for(i=0; i<nuse*nuse; i++) cov[i] = d2f[i];
   i = invert(nuse, cov, &chip);
   Mvec(nuse, cov, df, eigdpar);	/* -Offset from bottom of parabola */
#ifdef EIGBLAB
   printf("dpar in eigenvector basis: \n");
   for(i=0; i<nuse; i++) printf(" %16.8Le", eigdpar[i]);
   printf("\n");
#endif
   Mtrans(nuse, eigvec, eigvec);
   Mvec(nuse, eigvec, eigdpar, dpar);	/* -Offset in coord basis */
   Mtrans(nuse, eigvec, eigvec);
   for(j=0; j<nuse; j++) par[on[j]] -= dpar[j];
#ifdef EIGBLAB
   printf("dpar in coord basis: \n");
   func(npar, par, usepar, &chip, 0, NULL, NULL);
   for(i=0; i<nuse; i++) printf(" %16.8Le", dpar[i]);
   printf("   %16.8Le\n", chip-chi0);
#endif

/* Rotate d2f back to coordinates from eigenvector basis */
   Mmult(nuse, eigvec, d2f, cov);
   Mtrans(nuse, eigvec, eigvec);
   Mmult(nuse, cov, eigvec, d2f);
   Mtrans(nuse, eigvec, eigvec);

#ifdef EIGBLAB
   printf("d2f/2 in coord basis: \n");
   for(j=0; j<nuse; j++) {
      for(i=0; i<=j; i++) printf(" %16.8Le", d2f[i+j*nuse]/2);
      printf("\n");
   }
#endif

/* Recalculate the covariance matrix */
   for(i=0; i<nuse*nuse; i++) cov[i] = d2f[i] / 2;
   i = invert(nuse, cov, &chip);

#ifdef EIGBLAB
   printf("cov in coord basis: \n");
   for(j=0; j<nuse; j++) {
      for(i=0; i<=j; i++) printf(" %16.8Le", cov[i+j*nuse]);
      printf("\n");
   }
#endif

   for(j=0; j<nuse; j++) cov[j+j*nuse] = cov[j+j*nuse] > 0 ? 
			    SQRT(cov[j+j*nuse]) : -SQRT(-cov[j+j*nuse]);
   for(j=1; j<nuse; j++) {
      for(i=0; i<j; i++) {
	 cov[i+j*nuse] = cov[j+i*nuse] = 
	    cov[i+j*nuse] / (cov[i+i*nuse]*cov[j+j*nuse]);
      }
   }

#ifdef EIGBLAB
   printf("rcov in coord basis: \n");
   for(j=0; j<nuse; j++) {
      for(i=0; i<j; i++) printf(" %16.8Lf", cov[i+j*nuse]);
      printf(" %16.8Le\n", cov[j+j*nuse]);
   }
#endif

   free(dpar);
   free(df);

   return;
}

#if 0		/* Test code, to be ditched at some point */
/* Tune up numerical covariance matrix using pairwise eigenvectors */
void coveigen(int npar, REAL *par, int *usepar, MINIFUNC func, REAL *cov)
{
   int i, j, di, dj, nuse, on[MAXPAR];
   REAL chi, chi0, e1, e2, a, d, b, tg1, cth, sth;
   REAL A, B, D, ds1, dt1, ds2, dt2, d1, d2;
   REAL dfds, dfdt, d2fds, d2fdt, d2fdsdt, df[9];
   REAL DFVAR;
   REAL *TMPCURV, *TMPCOV;
      
   DFVAR = 0.01;
//   DFVAR = 1e-4;

/* Which parameters are actively varied? */
   for(i=0, nuse=0; i<npar; i++) {
      if(usepar[i]) {
	 on[nuse] = i;
	 nuse++;
      }
   }

   TMPCURV = (REAL *)calloc(nuse*nuse, sizeof(REAL));
   TMPCOV = (REAL *)calloc(nuse*nuse, sizeof(REAL));

/* Undo the "reduced covariance", invert to restore C^-1 = d2f/2 */
   for(j=1; j<nuse; j++) {
      for(i=0; i<j; i++) cov[i+nuse*j] = cov[j+nuse*i] = 
			    cov[i+nuse*j] * cov[i+nuse*i]*cov[j+nuse*j];
   }
   for(j=0; j<nuse; j++) {
      cov[j+nuse*j] = cov[j+nuse*j] > 0 ? cov[j+nuse*j]*cov[j+nuse*j] :
	                                 -cov[j+nuse*j]*cov[j+nuse*j];
   }
   invert(nuse, cov, &chi);

/* Center point for par */
   func(npar, par, usepar, &chi0, 0, NULL, NULL);

   printf("\n\n coveigen:\n");
   printf("%13.6Lf %13.6Lf %13.6Lf  %13.9Lf %13.9Lf %13.9Lf  %13.9Lf\n", 
	  par[on[0]], par[on[1]], par[on[2]], 
	  par[on[3]], par[on[4]], par[on[5]], chi0);
   printf("\n");

/* Initialize on-diagonal entries */
//   for(j=0; j<nuse; j++) cov[j+j*nuse] = 0;
   for(j=0; j<nuse; j++) TMPCURV[j+j*nuse] = 0;

/* Loop over all pairs of coordinates */
   for(j=1; j<nuse; j++) {
      for(i=0; i<j; i++) {

/* Hessian entries (divided by 2, ala mini convention) */
	 a = cov[i+i*nuse];
	 b = cov[i+j*nuse];
	 d = cov[j+j*nuse];
/* Eigenvalues */
	 e1 = (a+d)/2 - SQRT((a-d)*(a-d)/4+b*b);
	 e2 = (a+d)/2 + SQRT((a-d)*(a-d)/4+b*b);

/* Eigenvectors (slope tangents) */
	 tg1 = b / (e1-d);
	 printf("%d %d  e1= %12.4Le e2= %12.4Le  %12.4Le  slope= %12.4Lf\n", 
		i, j, e1, e2, e2/e1, tg1);

/* Sine and cosine */
	 cth = SQRT(1/(1+tg1*tg1));
	 sth = tg1 * cth;

/* Rotate to eigenvector basis: (c s / -s c) (a b / b d) (c -s / s c) */
/* Small eigenvalue along s */
	 A = a*cth*cth + 2*b*sth*cth + d*sth*sth;
	 D = a*sth*sth - 2*b*sth*cth + d*cth*cth;
	 B = (d-a)*sth*cth + b*(cth*cth-sth*sth);
//   printf("%12.4Le %12.4Le\n", A, B);
//   printf("%12.4Le %12.4Le\n", B, D);

/* Offsets to increase chi by DFVAR from the minimum */
	 d1 = SQRT(DFVAR / A);
	 d2 = SQRT(DFVAR / D);

/* Offsets along eigenvectors according to curvature */
	 ds1 =  d1*cth;
	 dt1 =  d1*sth;
	 ds2 = -d2*sth;
	 dt2 =  d2*cth;

/* Compute differences at 5 points (plus center) */
	 for(dj=-1; dj<=1; dj++) {
	    for(di=-1; di<=1; di++) {
	       if(di+1+(dj+1)*3 == 0 || di+1+(dj+1)*3 == 2 ||
		  di+1+(dj+1)*3 == 4 || di+1+(dj+1)*3 == 6) continue;
	       par[on[i]] += di*ds1 + dj*ds2;
	       par[on[j]] += di*dt1 + dj*dt2;
	       func(npar, par, usepar, &chi, 0, NULL, NULL);
//	       printf("%13.6Lf %13.6Lf %13.6Lf  %13.9Lf %13.9Lf %13.9Lf  %3d %3d %13.9Lf\n", 
//		      par[on[0]], par[on[1]], par[on[2]], 
//		      par[on[3]], par[on[4]], par[on[5]], di, dj, chi-chi0);
	       df[di+1+(dj+1)*3] = chi - chi0;
	       par[on[i]] -= di*ds1 + dj*ds2;
	       par[on[j]] -= di*dt1 + dj*dt2;
	    }
	 }

	 dfds = (df[5]-df[3]) / 2 / d1;
	 dfdt = (df[7]-df[1]) / 2 / d2;
	 d2fds = (df[5]+df[3]) / (d1*d1);
	 d2fdt = (df[7]+df[3]) / (d2*d2);
	 d2fdsdt = (df[8] - dfds*d1 - dfdt*d2 -
		    0.5*d2fds*d1*d1 - 0.5*d2fdt*d2*d2) / (d1*d2);
//	 printf("%14.6Le %14.6Le  %14.6Le %14.6Le %14.6Le\n", dfds, dfdt,
//		d2fds, d2fdt, d2fdsdt);

/* Rotate back */
	 A = d2fds/2*cth*cth + 2*d2fdsdt/2*sth*cth + d2fdt/2*sth*sth;
	 D = d2fds/2*sth*sth + 2*d2fdsdt/2*sth*cth + d2fdt/2*cth*cth;
	 B = -(d2fdt/2-d2fds/2)*sth*cth + d2fdsdt/2*(cth*cth-sth*sth);
//	 printf("%12.4Le %12.4Le\n", A, B);
//	 printf("%12.4Le %12.4Le\n", B, D);

//	 cov[i+nuse*i] += A;
	 cov[i+nuse*j] = cov[j+nuse*i] = B;
//	 cov[j+nuse*j] += D;

	 TMPCURV[i+nuse*i] += A;
	 TMPCURV[j+nuse*j] += D;

	 TMPCURV[j+nuse*i] += A;
	 TMPCURV[i+nuse*j] += D;

      }
   }

/* Average on-diagonal curvatures */
   for(j=0; j<nuse; j++) TMPCURV[j+j*nuse] /= nuse - 1;
   for(j=0; j<nuse; j++) cov[j+j*nuse] = TMPCURV[j+j*nuse];
   for(j=0; j<nuse; j++) TMPCURV[j+j*nuse] = cov[j+j*nuse];

/* Recalculate the covariance matrix */
   for(i=0; i<nuse*nuse; i++) TMPCOV[i] = cov[i] / 2;
   i = invert(nuse, TMPCOV, &chi);
   for(j=0; j<nuse; j++) TMPCOV[j+j*nuse] = TMPCOV[j+j*nuse] > 0 ? 
			    SQRT(TMPCOV[j+j*nuse]) : -SQRT(-TMPCOV[j+j*nuse]);
   for(j=1; j<nuse; j++) {
      for(i=0; i<j; i++) {
	 TMPCOV[i+j*nuse] = TMPCOV[j+i*nuse] = 
	    TMPCOV[i+j*nuse] / (TMPCOV[i+i*nuse]*TMPCOV[j+j*nuse]);
      }
   }

/* Tell us about everything: d2f, values for curvature, covariance */
   printf("\n curv:\n");
   for(j=0; j<nuse; j++) {
      for(i=0; i<nuse; i++) printf(" %12.4Le", TMPCURV[i+j*nuse]);
      printf("\n");
   }

   printf("\n d2f:\n");
   for(j=0; j<nuse; j++) {
      for(i=0; i<=j; i++) printf(" %12.4Le", cov[i+j*nuse]);
      printf("\n");
   }

/* Recalculate the covariance matrix */
   for(i=0; i<nuse*nuse; i++) TMPCOV[i] = cov[i] / 2;
   i = invert(nuse, TMPCOV, &chi);

   printf("\n inv d2f:\n");
   for(j=0; j<nuse; j++) {
      for(i=0; i<=j; i++) printf(" %12.4Le", TMPCOV[i+j*nuse]);
      printf("\n");
   }

   for(j=0; j<nuse; j++) TMPCOV[j+j*nuse] = TMPCOV[j+j*nuse] > 0 ? 
			    SQRT(TMPCOV[j+j*nuse]) : -SQRT(-TMPCOV[j+j*nuse]);
   for(j=1; j<nuse; j++) {
      for(i=0; i<j; i++) {
	 TMPCOV[i+j*nuse] = TMPCOV[j+i*nuse] = 
	    TMPCOV[i+j*nuse] / (TMPCOV[i+i*nuse]*TMPCOV[j+j*nuse]);
      }
   }

   printf("\n cov:\n");
   for(j=0; j<nuse; j++) {
      for(i=0; i<=j; i++) printf(" %12.4Le", TMPCOV[i+j*nuse]);
      printf("\n");
   }

   return;
}
#endif
