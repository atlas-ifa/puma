/* A few sorting routines */
/* Subversion r.13388 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <string>
#include <math.h>
#include "tsort.h"

#define QSORT_NMIN 100  /* insort instead of qsort if fewer than this  */
#define MAXSTACK 256	/* internal storage for qsort */
#define NSTOP 10	/* out of order limit for qsort swapping */

#define MEDEC_BIN        10     /* quantile decimation factor */
#define MEDEC_NOSIFT    200	/* fewer: don't decimate at all */
#define MEDEC_ONEPASS  1000	/* fewer: decimate once */
#define MEDEC_TWOPASS 20000	/* fewer: decimate twice */
#define INSHEAP 200             /* transition from insertion to heap median */

// #define VERIFY_SORT	/* Verify every sort */

// #define USE_MEMMOVE  /* you feeling lucky, punk? */

/* Pseudo-random number algorithm of Marsaglia: rand = xorshift32(&rand) */
uint32_t xorshift32(uint32_t *state)
{
   uint32_t x = *state;
   x ^= x << 13;
   x ^= x >> 17;
   x ^= x << 5;
   return x;
}

/* Return the number of sorted x such that x[k] < key */
template <typename DATA>
int verify_sort(int n, DATA *x)
{
   register int k, order;
   for(k=0, order=0; k<n-1; k++) order += x[k] <= x[k+1];
   if(order < n-1) fprintf(stderr, "Verify error: %d  ", order);
   return(order==n-1);
}

/* Return the number of sorted x such that x[k] < key */
int verify_sort_str(int n, char **x)
{
   register int k, order;
   for(k=0, order=0; k<n-1; k++) order += strcmp(x[k],x[k+1]) <= 0;
   if(order < n-1) fprintf(stderr, "Verify error: %d  ", order);
   return(order==n-1);
}

/* Return the number of sorted x such that x[k] < key */
template <typename DATA>
int linlt(int n, DATA *x, DATA key)
{
   register int k;
   for(k=0; k<n; k++) if(x[k] >= key) return(k);
   return(n);
}

/* Return the number of sorted x such that x[k] < key */
template <typename DATA>
int binlt(int n, DATA *x, DATA key)
{
   register int k0=0, k1=n-1, k;
   if(key <= x[k0] ) return(0);
   if(key >  x[k1] ) return(n);
   while(k1 > k0+1) {  // chop, maintaining x[k0] < key and x[k1] >= key
      k = (k1 + k0) / 2;
      if(x[k] < key) k0 = k;
      else           k1 = k;
   }
   return(k1);
}

/* Insertion sort, possibly carry an auxiliary array 020130 John Tonry */
/* Binary chop for randomly distributed data */
template <typename DATA, typename AUX>
void insort2(int n, DATA *x, AUX *idx)
{
   register int j, k;
   register DATA tmp;
   register AUX itmp;
#ifndef USE_MEMMOVE
   register int i;
#endif

   if(n < 2) return;
   if(idx == NULL) {
      for(j=n-2; j>=0; j--) {
	 k = j + binlt<DATA>(n-j-1, x+j+1, x[j]);

#ifdef USE_MEMMOVE
	 if(k != j) {
	    tmp = x[j];
	    memmove(x+j, x+j+1, (k-j)*sizeof(DATA));
	    x[k] = tmp;
	 }
#else
	tmp = x[j];
	for(i=j+1; i<=k; i++) x[i-1] = x[i];
	x[k] = tmp;
#endif

      }
   } else {
      for(j=n-2; j>=0; j--) {
	 k = j + binlt<DATA>(n-j-1, x+j+1, x[j]);
	 if(k != j) {
#ifdef USE_MEMMOVE
	    tmp = x[j];
	    memmove(x+j, x+j+1, (k-j)*sizeof(DATA));
	    x[k] = tmp;
	    itmp = idx[j];
	    memmove(idx+j, idx+j+1, (k-j)*sizeof(AUX));
	    idx[k] = itmp;
#else
            tmp = x[j];
            itmp = idx[j];
            for(i=j+1; i<=k; i++) {
               x[i-1] = x[i];
               idx[i-1] = idx[i];
            }
            x[k] = tmp;
            idx[k] = itmp;
#endif
	 }
      }
   }
#ifdef VERIFY_SORT
   j = verify_sort<DATA>(n, x);
#endif
   return;
}

/* Insertion sort, possibly carry an auxiliary array */
/* Linear search for arrays that are nearly sorted */
template <typename DATA, typename AUX>
void linsort2(int n, DATA *x, AUX *idx)
{
   register int i, j, k;
   register DATA tmp;
   register AUX itmp;

   if(n < 2) return;
   if(idx == NULL) {
      for(j=n-2; j>=0; j--) {
	 for(i=j+1, k=j; i<n; i++) {
	    if(x[j] <= x[i]) break;
	    k = i;
	 }
	 if(k != j) {
#ifdef USE_MEMMOVE
	    tmp = x[j];
	    memmove(x+j, x+j+1, (k-j)*sizeof(DATA));
	    x[k] = tmp;
#else
            tmp = x[j];
            for(i=j+1; i<=k; i++) x[i-1] = x[i];
            x[k] = tmp;
#endif
	 }
      }
   } else {
      for(j=n-2; j>=0; j--) {
	 for(i=j+1, k=j; i<n; i++) {
	    if(x[j] <= x[i]) break;
	    k = i;
	 }
	 if(k != j) {
#ifdef USE_MEMMOVE
	    tmp = x[j];
	    memmove(x+j, x+j+1, (k-j)*sizeof(DATA));
	    x[k] = tmp;
	    itmp = idx[j];
	    memmove(idx+j, idx+j+1, (k-j)*sizeof(AUX));
	    idx[k] = itmp;
#else
            tmp = x[j];
            itmp = idx[j];
            for(i=j+1; i<=k; i++) {
               x[i-1] = x[i];
               idx[i-1] = idx[i];
            }
            x[k] = tmp;
            idx[k] = itmp;
#endif
	 }
      }
   }
   return;
}

/* Quicksort, possibly carry an auxiliary array: 1980 (John Tonry) */
template <typename DATA, typename AUX>
void qsort2(int n, DATA *x, AUX *idx)
{
   register DATA key, temp;
   register AUX itemp;
   int l, r, m, lstack[MAXSTACK], rstack[MAXSTACK], randomness=1, sp;
   uint32_t ran=3141592653;
   register int i, j;

   if(n < QSORT_NMIN) goto INSORT;	// Insertion sorting when n is small

   for(i=randomness=0; i<n-1; i++) randomness += x[i+1] > x[i];
   if(randomness > n/2) randomness = n - randomness;
   if(randomness < n/2 - 5*sqrt((double)n)) randomness = 0;

   sp = 0;
   lstack[sp] = 0;
   rstack[sp] = n-1;

   while(sp >= 0) {		// Sort a subrecord off the stack
      l = lstack[sp];
      r = rstack[sp--];

      if(!randomness) {	// Set KEY = random point between X(L), X(R)
	 ran = xorshift32(&ran);
	 m = l + ran / 4294967296.0 * (r-l);
	 key = x[m];

      } else {			// Set KEY = median of X(L), X(M), X(R)
         m = (l + r) / 2;
	 if((x[m]>x[l]) ^ (x[r]>x[m])) {
	    if((x[m]>x[l]) ^ (x[l]>x[r])) key = x[r];
	    else                          key = x[l];
	 } else {
	    key = x[m];
	 }
      }

      i = l;
      j = r;
      if(idx == NULL) {
	 while(1) {
	    while(x[i] < key) i++;	// Find a big record on the left
	    while(x[j] > key) j--;	// Find a small record on the right
	    if(i >= j) break;
	    temp = x[i];	    x[i++] = x[j];	    x[j--] = temp;
	 }
      } else {
	 while(1) {
	    while(x[i] < key) i++;	// Find a big record on the left
	    while(x[j] > key) j--;	// Find a small record on the right
	    if(i >= j) break;
	    temp = x[i];	    x[i] = x[j];	    x[j] = temp;
	    itemp = idx[i];	    idx[i++] = idx[j];	    idx[j--] = itemp;
	 }
      }

/* Subfile is partitioned into two halves, left .le. right */
      if(j-l+1 > NSTOP) {	// Push the two halves on the stack
	 if(sp >= MAXSTACK) break;	// Stack overflow, default to insertion
	 lstack[++sp] = l;
	 rstack[sp] = j;
      }
      if(r-j > NSTOP) {
	 if(sp >= MAXSTACK) break;	// Stack overflow, default to insertion
	 lstack[++sp] = j+1;
	 rstack[sp] = r;
      }
   }

INSORT:					// Insertion sorting to finish up
   linsort2<DATA, AUX>(n, x, idx);
#ifdef VERIFY_SORT
   j = verify_sort<DATA>(n, x);
#endif
   return;
}


/* decimate() moves central points which straddle goal to start of array */
/* There are nnew points left in the array that must be sorted, and the */
/* goal parameter is updated to where it lies in the subarray */
/* Typical usage:
 *
 *    q = 0.5;
 *    decimate(n, x, 10, &q, &nnew);
 *    qsort2(nnew, x, NULL);
 *    k = q*(nnew-1) - 0.49999999999;
 *    frac = q*(nnew-1) - k;
 *    quantile = (1-frac)*x[k] + frac*x[k+1];
 *
 * Note that decimate() can be called multiple times and n can be updated:
 *
 *    n = nx;
 *    q = 0.5;
 *    decimate(n, x, 10, &q, &n);
 *    decimate(n, x, 10, &q, &n);
 *    qsort2(n, x, NULL);
 *    k = q*(n-1) - 0.49999999999;
 *    frac = q*(n-1) - k;
 *    quantile = (1-frac)*x[k] + frac*x[k+1];
 */
template <typename DATA>
void decimate(int n, DATA *x, int nbin, double *goal, int *nnew)
{
   int i, j, k, k0, *sift_ngelt, nbelow;
   DATA xtmp, *sift_posts;

   if(n < 4*nbin) {
      *nnew = n;
      return;
   }

/* Both x[k0] and x[k0+1] must exist in the decimated set */
   k0 = (*goal)*(n-1) - 0.49999999999;
   if(k0 < 0) k0 = 0;

/* Allocate some buffer space */
   sift_posts = (DATA *)calloc(nbin+1, sizeof(DATA));
   sift_ngelt = (int *)calloc(nbin, sizeof(int));	// Note init to 0

/* Sample some points */
   for(i=1; i<nbin; i++) sift_posts[i] = x[(i*(n-1))/nbin];
/* and sort them */
   insort2<DATA, void*>(nbin-1, sift_posts+1, NULL);

/* Here's what we're calculating
 *    sift_posts[0]           = minimum value
 *    sift_posts[1:nbin-1]    = sorted division points among the values
 *    sift_posts[nbin]        = maximum value
 *    
 *    sift_ngelt[0]    = number of x: ( x <= sift_posts[1] )
 *    sift_ngelt[j]    = number of x: ( sift_posts[j] < x <= sift_posts[j+1] )
 *       i.e. x[i] == sift_posts[j=1:nbin-1] increments sift_ngelt[j-1]
 */

   sift_posts[0] = sift_posts[1];
   sift_posts[nbin] = sift_posts[nbin-1];
/* Count up all the points and find the min and max */
   for(i=0; i<n; i++) {
      j = binlt<DATA>(nbin-1, sift_posts+1, x[i]);	// NB: base addr + 1
      sift_ngelt[j] += 1;
      if(j == 0 && x[i] < sift_posts[0]) sift_posts[0] = x[i];
      if(j == nbin-1 && x[i] > sift_posts[nbin]) sift_posts[nbin] = x[i];
   }

/* Figure out the points which can be ignored for the median */
/*   frac = (*goal)*(n-1) - k0; */
/*   quantile = (1-frac)*x[k0] + frac*x[k0+1]; */
/* The goal lies between sift_posts[j] and sift_posts[k] */
   for(j=0, nbelow=0; j<nbin; j++) {
      if(nbelow + sift_ngelt[j] > k0) break;
      nbelow += sift_ngelt[j];
   }
   k = (nbelow + sift_ngelt[j] > k0+1) ? j+1 : j+2;

/* Swap the points between sift_posts[j] and sift_posts[k] to the beginning */
   if(j == 0) {
      for(i=0, *nnew=0; i<n; i++) {
	 if(x[i] <= sift_posts[k]) {
	    xtmp = x[i];
	    x[i] = x[*nnew];
	    x[*nnew] = xtmp;
	    *nnew += 1;
	 }
      }
   } else if(k == nbin) {
      for(i=0, *nnew=0; i<n; i++) {
	 if(sift_posts[j] < x[i]) {
	    xtmp = x[i];
	    x[i] = x[*nnew];
	    x[*nnew] = xtmp;
	    *nnew += 1;
	 }
      }
   } else {
      for(i=0, *nnew=0; i<n; i++) {
	 if(sift_posts[j] < x[i] && x[i] <= sift_posts[k]) {
	    xtmp = x[i];
	    x[i] = x[*nnew];
	    x[*nnew] = xtmp;
	    *nnew += 1;
	 }
      }
   }

/* Update the goal for this diminished set */
   *goal = ((*goal)*(n-1) - nbelow) / (*nnew-1);

   return;
}

/* 090315 John Tonry */
template <typename DATA>
void wdecimate(int n, DATA *x, double *wgt, int nbin, double *goal, int *nnew)
{
   int i, j, k, loend, hiend;
   DATA xtmp, lo, hi, xhi;
   double tmp, below, wtot;
   DATA *sift_posts;
   double *sift_gelt;

   if(n < nbin+1) {
      *nnew = n;
      below = 0;
      return;
   }

/* Allocate some buffer space */
   sift_posts = (DATA *)calloc(nbin+1, sizeof(DATA));
   sift_gelt = (double *)calloc(nbin, sizeof(double));	// Note init to 0

/* Sample some points */
   for(i=1; i<nbin; i++) sift_posts[i] = x[(i*(n-1))/nbin];
/* and sort them */
   insort2<DATA, void*>(nbin-1, sift_posts+1, NULL);

/* Here's what we're calculating
 *    sift_posts[0]           = minimum value
 *    sift_posts[1:nbin-1]    = sampled division points among the values
 *    sift_posts[nbin]        = maximum value
 *    
 *    sift_gelt[0]    = sum weight: ( x <= sift_posts[1] )
 *    sift_gelt[j]    = sum weight: ( sift_posts[j] < x <= sift_posts[j+1] )
 *       i.e. x[i] == sift_posts[j=1:nbin-1] augments sift_gelt[j-1]
 */

   sift_posts[0] = sift_posts[1];
   sift_posts[nbin] = sift_posts[nbin-1];
/* Count up all the points and find the min and max */
   for(i=0, wtot=0.0; i<n; i++) {
      wtot += wgt[i];
      j = binlt<DATA>(nbin-1, sift_posts+1, x[i]);	// NB: base addr + 1
      sift_gelt[j] += wgt[i];
      if(j == 0 && x[i] < sift_posts[0]) sift_posts[0] = x[i];
      if(j == nbin-1 && x[i] >= sift_posts[nbin]) sift_posts[nbin] = x[i];
   }

/* Figure out the points which can be ignored */
/* The points between sift_posts[j] and sift_posts[j+1] contain the goal */
   for(j=0, below=0.0; j<nbin; j++) {
      if(below + sift_gelt[j] >= (*goal)*wtot) break;
      below += sift_gelt[j];
   }

/* Swap the points between sift_posts[j] and sift_posts[k] to the beginning */
/* Test whether the goal is achieved at an endpoint */
   loend = hiend = 0;
   lo = x[0];
   hi = x[n-1];
   if(j == 0) {
      for(i=0, *nnew=0; i<n; i++) {
	 if(x[i] <= sift_posts[j+1]) {
	    if(below + wgt[i] > (*goal)*wtot) loend = 1;
	    if(below + sift_gelt[j] - wgt[i] < (*goal)*wtot) hiend = 1;
	    if(*nnew == 0) {
	       lo = hi = x[i];
	    } else {
	       if(x[i] < lo) lo = x[i];
	       if(x[i] > hi) hi = x[i];
	    }
	    xtmp = x[i];    x[i] = x[*nnew];     x[*nnew] = xtmp;
	    tmp = wgt[i]; wgt[i] = wgt[*nnew]; wgt[*nnew] = tmp;
	    *nnew += 1;
	 }
      }
   } else if(j == nbin-1) {
      for(i=0, *nnew=0; i<n; i++) {
	 if(sift_posts[j] < x[i]) {
	    if(below + wgt[i] > (*goal)*wtot) loend = 1;
	    if(below + sift_gelt[j] - wgt[i] < (*goal)*wtot) hiend = 1;
	    if(*nnew == 0) {
	       lo = hi = x[i];
	    } else {
	       if(x[i] < lo) lo = x[i];
	       if(x[i] > hi) hi = x[i];
	    }
	    xtmp = x[i];    x[i] = x[*nnew];     x[*nnew] = xtmp;
	    tmp = wgt[i]; wgt[i] = wgt[*nnew]; wgt[*nnew] = tmp;
	    *nnew += 1;
	 }
      }
   } else {
      for(i=0, *nnew=0; i<n; i++) {
	 if(sift_posts[j] < x[i] && x[i] <= sift_posts[j+1]) {
	    if(below + wgt[i] > (*goal)*wtot) loend = 1;
	    if(below + sift_gelt[j] - wgt[i] < (*goal)*wtot) hiend = 1;
	    if(*nnew == 0) {
	       lo = hi = x[i];
	    } else {
	       if(x[i] < lo) lo = x[i];
	       if(x[i] > hi) hi = x[i];
	    }
	    xtmp = x[i];    x[i] = x[*nnew];     x[*nnew] = xtmp;
	    tmp = wgt[i]; wgt[i] = wgt[*nnew]; wgt[*nnew] = tmp;
	    *nnew += 1;
	 }
      }
   }

/* Augment the new set if loend or hiend is set */
   if(loend) {		// goal is met by the smallest entry
      k = 0;
      xhi = sift_posts[0];
      for(i=*nnew; i<n; i++) {
	 if(x[i] <= lo && x[i] >= xhi) {	// biggest but smaller than lo
	    xhi = x[i];
	    k = i;
	 }
      }
      if(k > 0) {
	 xtmp = x[k];    x[k] = x[*nnew];      x[*nnew] = xtmp;
	 tmp = wgt[k]; wgt[k] = wgt[*nnew];  wgt[*nnew] = tmp;
	 *nnew += 1;
	 below -= wgt[k];
	 sift_gelt[j] += wgt[k];
      }
   }

   if(hiend) {		// goal is met by the largest entry
      k = 0;
      xhi = sift_posts[nbin];
      for(i=*nnew; i<n; i++) {
	 if(x[i] >= hi && x[i] <= xhi) {	// smallest but bigger than hi
	    xhi = x[i];
	    k = i;
	 }
      }
      if(k > 0) {
	 xtmp = x[k];     x[k] = x[*nnew];      x[*nnew] = xtmp;
	 tmp = wgt[k];  wgt[k] = wgt[*nnew];  wgt[*nnew] = tmp;
	 *nnew += 1;
	 sift_gelt[j] += wgt[k];
      }
   }

/* Update goal for this new set */
   *goal = (*goal*wtot - below) / sift_gelt[j];

   free(sift_posts);
   free(sift_gelt);

   return;
}


/* Return weighted quantile */
template <typename DATA>
DATA wquant(int n, DATA *x, double *wgt, double goal)
//   int n,       = number of data points
//   DATA *x,     = values of data points: will be sorted
//   double *wgt, = weights of data points wgt[idx[i]] corresponds to x[i]
//   double goal  = Return this quantile: 0.5 for median
{
   int i, *idx;
   DATA xm0, xm1;
   double wtot, wprev, wthis, wm0, wm1;

/* Degenerate cases */
   if(n <= 0) {
      return((DATA)0);
   } else if(n == 1) {
      return(x[0]);
   }

/* sort x and index arrays */
   idx = (int *)calloc(n, sizeof(int));
   for(i=0, wtot=0.0; i<n; i++) {
      idx[i] = i;
      wtot += wgt[i];
   }
   qsort2<DATA, int>(n, x, idx);

/* All zero weights?  Naughty, return something plausible. */
   if(wtot <= 0) {
      free(idx);
      return(x[n/2]);
   }

/* Find the point which straddles the desired goal */
   wprev = 0.0;
   xm0 = x[0];
   xm1 = x[n-1];
   wm0 = 0.0;
   wm1 = wtot;
   for(i=0; i<n; i++) {
      wthis = wgt[idx[i]];
      if(wprev < goal*wtot && wprev+wthis >= goal*wtot) {
	 if(i > 0)   xm0 = 0.5 * ( x[i-1] +   x[i]);
	 else        xm0 = 0.5 * ( 3*x[i] -   x[i+1]);
	 if(i < n-1) xm1 = 0.5 * (   x[i] +   x[i+1]);
	 else        xm1 = 0.5 * (-x[i-1] + 3*x[i]);
	 wm0 = wprev / wtot;
	 wm1 = (wprev + wthis) / wtot;
	 break;
      }
      wprev += wthis;
   }
   free(idx);
   return((DATA)(xm0 + (xm1-xm0) * (goal-wm0)/(wm1-wm0)));
}





/* Find a quantile q of data x, possibly using weights */
template <typename DATA>
DATA quantile(int n, DATA *x, double *wgt, double q)
{
   double frac;
   int goal, goal1;
   DATA quantile;

   quantile = 0;
   if(n <= 0) return((DATA)0);
   if(n == 1) return(x[0]);

   if(q > 1) q = 1;
   if(q < 0) q = 0;

   if(wgt == NULL) {
      if(n >= MEDEC_TWOPASS) decimate<DATA>(n, x, MEDEC_BIN, &q, &n);
      if(n >= MEDEC_ONEPASS) decimate<DATA>(n, x, MEDEC_BIN, &q, &n);
      if(n >= MEDEC_NOSIFT)  decimate<DATA>(n, x, MEDEC_BIN, &q, &n);

      qsort2<DATA, void*>(n, x, NULL);

      goal = q*(n-1) + 0.49999999999;
      goal1 = goal<n-1 ? goal+1 : goal;
      frac = q*(n-1) - goal;

      quantile = (1-frac)*x[goal] + frac*x[goal1];

   } else {
      if(n >= MEDEC_TWOPASS) wdecimate<DATA>(n, x, wgt, MEDEC_BIN, &q, &n);
      if(n >= MEDEC_ONEPASS) wdecimate<DATA>(n, x, wgt, MEDEC_BIN, &q, &n);
      if(n >= MEDEC_NOSIFT)  wdecimate<DATA>(n, x, wgt, MEDEC_BIN, &q, &n);

      quantile = wquant<DATA>(n, x, wgt, q);
   }

   return(quantile);
}


/* Quicksort string pointers, possibly carry an auxiliary array */
template <typename AUX>
void qsort2_str(int n, char **x, AUX *idx)
{
   register char *key, *temp;
   register AUX itemp;
   int l, r, m, lstack[MAXSTACK], rstack[MAXSTACK], randomness=1, sp;
   uint32_t ran=3141592653;
   register int i, j, k;

   if(n < QSORT_NMIN) goto INSORT;	// Insertion sorting when n is small

   for(i=randomness=0; i<n-1; i++) randomness += strcmp(x[i+1],x[i]) > 0;
   if(randomness > n/2) randomness = n - randomness;
   if(randomness < n/2 - 5*sqrt((double)n)) randomness = 0;

   sp = 0;
   lstack[sp] = 0;
   rstack[sp] = n-1;

   while(sp >= 0) {		// Sort a subrecord off the stack
      l = lstack[sp];
      r = rstack[sp--];

      if(!randomness) {	// Set KEY = random point between X(L), X(R)
	 ran = xorshift32(&ran);
	 m = l + ran / 4294967296.0 * (r-l);
	 key = x[m];

      } else {			// Set KEY = median of X(L), X(M), X(R)
         m = (l + r) / 2;
	 if((strcmp(x[m],x[l])>0) ^ (strcmp(x[r],x[m])>0)) {
	    if((strcmp(x[m],x[l])>0) ^ (strcmp(x[l],x[r])>0)) key = x[r];
	    else                                              key = x[l];
	 } else {
	    key = x[m];
	 }
      }

      i = l;
      j = r;
      if(idx == NULL) {
	 while(1) {
	    while(strcmp(x[i],key) < 0) i++; // Find a big record on the left
	    while(strcmp(x[j],key) > 0) j--; // Find a small record on the right
	    if(i >= j) break;
	    temp = x[i];	    x[i++] = x[j];	    x[j--] = temp;
	 }
      } else {
	 while(1) {
	    while(strcmp(x[i],key) < 0) i++; // Find a big record on the left
	    while(strcmp(x[j],key) > 0) j--; // Find a small record on the right
	    if(i >= j) break;
	    temp = x[i];	    x[i] = x[j];	    x[j] = temp;
	    itemp = idx[i];	    idx[i++] = idx[j];	    idx[j--] = itemp;
	 }
      }

/* Subfile is partitioned into two halves, left .le. right */
      if(j-l+1 > NSTOP) {	// Push the two halves on the stack
	 if(sp >= MAXSTACK) break;	// Stack overflow, default to insertion
	 lstack[++sp] = l;
	 rstack[sp] = j;
      }
      if(r-j > NSTOP) {
	 if(sp >= MAXSTACK) break;	// Stack overflow, default to insertion
	 lstack[++sp] = j+1;
	 rstack[sp] = r;
      }
   }

INSORT:					// Insertion sorting to finish up
   for(j=n-2; j>=0; j--) {
      for(i=j+1, k=j; i<n; i++) {
	 if(strcmp(x[j],x[i]) <= 0) break;
	 k = i;
      }
      if(k != j) {
#ifdef USE_MEMMOVE
	 temp = x[j];
	 memmove(x+j, x+j+1, (k-j)*sizeof(char *));
	 x[k] = temp;
	 if(idx != NULL) {
	    itemp = idx[j];
	    memmove(idx+j, idx+j+1, (k-j)*sizeof(AUX));
	    idx[k] = itemp;
	 }
#else
	temp = x[j];
	itemp = idx[j];
	for(i=j+1; i<=k; i++) {
	   x[i-1] = x[i];
	   idx[i-1] = idx[i];
	}
	x[k] = temp;
	idx[k] = itemp;
#endif
      }
   }
#ifdef VERIFY_SORT
   j = verify_sort_str(n, x);
#endif

   return;
}

////////////////////////
// Just a fast median //
////////////////////////

/* Network for median of 3 from http://www.cs.hut.fi/~cessu/selection */
#define M3_01(a,b,c) (a > b) ? M3_02(a,b,c) : M3_02(b,a,c)
#define M3_02(a,b,c) (a > c) ? M3_03(b,c) : a
#define M3_03(a,b) (a > b) ? a : b

/* Network for median of 4 (JT) */
#define M4_01(a,b,c,d) (a > b) ? M4_02(a,b,c,d) : M4_02(b,a,c,d) // a b c d
#define M4_02(a,b,c,d) (c > d) ? M4_03(a,c,b,d) : M4_03(a,d,b,c) // a>b    
#define M4_03(a,b,c,d) (a > b) ? M4_04(d,b,c) : M4_04(c,a,d)     // a>c b>d
#define M4_04(a,b,c) (a < c) ? (b+c)/2 : (a+b)/2                 // b>c

/* Network for median of 5 from http://www.cs.hut.fi/~cessu/selection */
#define M5_01(a,b,c,d,e) (a > b) ? M5_02(a,b,c,d,e) : M5_02(b,a,c,d,e) // a b c d e
#define M5_02(a,b,c,d,e) (c > d) ? M5_03(a,c,b,d,e) : M5_03(a,d,b,c,e) // a>b c d e
#define M5_03(a,b,c,d,e) (a > b) ? M5_04(b,d,c,e) : M5_04(a,c,d,e)     // a>c b>d e
#define M5_04(a,b,c,d) (c > d) ? M5_05(a,c,b,d) : M5_05(a,d,b,c)       // a>b      
#define M5_05(a,b,c,d) (a > b) ? M5_06(b,c) : M5_06(a,d)               // a>c b>d  
#define M5_06(a,b) (a > b) ? a : b

/* Network for median of 6 (JT) */
#define M6_01(a,b,c,d,e,f) (a > b) ? M6_02(a,b,c,d,e,f) : M6_02(b,a,c,d,e,f) // a b c d e f
#define M6_02(a,b,c,d,e,f) (c > d) ? M6_03(a,b,c,d,e,f) : M6_03(a,b,d,c,e,f) // a>b c d e f
#define M6_03(a,b,c,d,e,f) (a > c) ? M6_04(a,b,c,d,e,f) : M6_04(c,d,a,b,e,f) // a>b c>d e f
#define M6_04(a,b,c,d,e,f) (a > e) ? M6_05(c,d,b,e,f) : M6_06(a,b,c,d,f)     // a>b a>c>d e f
#define M6_05(a,b,c,d,e) (c > d) ? M6_07(a,b,c,d,e) : M6_07(a,b,d,c,f)       // a>b c d e H1
#define M6_06(a,b,c,d,e) (a > e) ? M6_08(c,d,b,e) : M6_10(a,b,c)             // a>b a>c>d e H1
#define M6_07(a,b,c,d,e) (a > c) ? M6_06(a,b,c,d,e) : M6_06(c,d,a,b,f)       // a>b c>d e H1
#define M6_08(a,b,c,d) (b > c) ? M6_10(a,b,d) : M6_09(a,c,d)                 // a>b c d H2
#define M6_09(a,b,c) (a > b) ? M6_10(a,b,c) : M6_10(b,a,c)                   // a b c H2,L1
#define M6_10(a,b,c) (b > c) ? (a+b)/2 : (a+c)/2                             // a>b c H2,L1

/* Network for median of 7 from http://www.cs.hut.fi/~cessu/selection */
#define M7_01(a,b,c,d,e,f,g) (a > b) ? M7_02(a,b,c,d,e,f,g) : M7_02(b,a,c,d,e,f,g) // a b c d e f g
#define M7_02(a,b,c,d,e,f,g) (c > d) ? M7_03(a,c,b,d,e,f,g) : M7_03(a,d,b,c,e,f,g) // a>b c d e f g
#define M7_03(a,b,c,d,e,f,g) (a > b) ? M7_04(a,b,d,c,e,f,g) : M7_04(b,a,c,d,e,f,g) // a>c,b>d
#define M7_04(a,b,c,d,e,f,g) (d > e) ? M7_05(b,d,c,e,f,g) : M7_06(a,b,e,c,d,f,g)   // a>b>c,a>d
#define M7_05(a,b,c,d,e,f) (a > b) ? M7_07(a,b,d,c,e,f) : M7_07(b,a,c,d,e,f)       // a>c,b>d
#define M7_06(a,b,c,d,e,f,g) (b > c) ? M7_07(b,c,e,d,f,g) : M7_08(a,c,b,d,e,f,g)   // a>b>d,a>e,c>e            
#define M7_07(a,b,c,d,e,f) (b > d) ? M7_16(a,b,c,d,e,f) : M7_17(a,d,b,e,f)         // a>b>c,a>d
#define M7_08(a,b,c,d,e,f,g) (f > g) ? M7_09(a,b,c,f,d,e,g) : M7_09(a,b,c,g,d,e,f) // a>c>d,a>e,b>c,b>d,b>e
#define M7_09(a,b,c,d,e,f,g) (c > g) ? M7_10(c,d,e,f) : M7_11(d,b,a,g,f)           // a>c>e,a>f,b>c,b>e,b>f,d>g
#define M7_10(a,b,c,d) (b > d) ? M7_14(a,b,c,d) : M7_14(a,d,c,b)                   // a>c
#define M7_11(a,b,c,d,e) (b > c) ? M7_12(c,a,e,d) : M7_12(b,a,e,d)                 // a>d,b>e,c>e
#define M7_12(a,b,c,d) (c > d) ? M7_13(c,b) : M7_13(d,a)                           // a>c,b>d
#define M7_13(a,b) (a > b) ? b : a                                                 // 
#define M7_14(a,b,c,d) (a > b) ? M7_15(b,c) : M7_15(a,d)                           // a>c,b>d
#define M7_15(a,b) (a > b) ? a : b                                                 // 
#define M7_16(a,b,c,d,e,f) (e > f) ? M7_20(a,b,e,d,c,f) : M7_20(a,b,f,d,c,e)       // a>b>c,a>d,b>d
#define M7_17(a,b,c,d,e) (a > d) ? M7_18(b,c,d,e) : M7_19(a,b,e)                   // a>b>c
#define M7_18(a,b,c,d) (c > d) ? M7_14(a,c,b,d) : M7_14(a,d,b,c)                   // a>b
#define M7_19(a,b,c)  (a > c) ? M7_15(b,c) : a                                     // a>b
#define M7_20(a,b,c,d,e,f) (b > c) ? M7_21(c,d,e) : M7_19(a,b,f)                   // a>b>d,a>e,b>e,c>f
#define M7_21(a,b,c) (a > b) ? M7_15(a,c) : M7_15(b,c)

/* Median of 3 numbers */
template <typename DATA>
DATA med3(DATA a, DATA b, DATA c)
{
   return(M3_01(a, b, c));
}

/* Median of 4 numbers (average middle 2) */
template <typename DATA>
DATA med4(DATA a, DATA b, DATA c, DATA d)
{
   return(M4_01(a, b, c, d));
}

/* Median of 5 numbers */
template <typename DATA>
DATA med5(DATA a, DATA b, DATA c, DATA d, DATA e)
{
   return(M5_01(a, b, c, d, e));
}

/* Median of 6 numbers (average middle 2) */
template <typename DATA>
DATA med6(DATA a, DATA b, DATA c, DATA d, DATA e, DATA f)
{
   return(M6_01(a, b, c, d, e, f));
}

/* Median of 7 numbers */
template <typename DATA>
DATA med7(DATA a, DATA b, DATA c, DATA d, DATA e, DATA f, DATA g)
{
   return(M7_01(a, b, c, d, e, f, g));
}

/* Median of 8 numbers (slower than insmed for g++ -O3) */
template <typename DATA>
DATA med8(DATA a, DATA b, DATA c, DATA d, DATA e, DATA f, DATA g, DATA h)
{
   register DATA tmp;
   if(a>e) {tmp=a; a=e; e=tmp;}
   if(b>f) {tmp=b; b=f; f=tmp;}
   if(c>g) {tmp=c; c=g; g=tmp;}
   if(d>h) {tmp=d; d=h; h=tmp;}
   if(a>c) {tmp=a; a=c; c=tmp;}
   if(b>d) {tmp=b; b=d; d=tmp;}
   if(a>b) {tmp=a; a=b; b=tmp;}
   if(e>g) {tmp=e; e=g; g=tmp;}
   if(f>h) {tmp=f; f=h; h=tmp;}
   if(g>h) {tmp=g; g=h; h=tmp;}
   return(med6(b, c, d, e, f, g));
}

/* Median of 9 numbers (slower than insmed for g++ -O3) */
template <typename DATA>
DATA med9(DATA a, DATA b, DATA c, DATA d, DATA e, DATA f, DATA g, DATA h, DATA i)
{
   register DATA tmp;
   if(a>e) {tmp=a; a=e; e=tmp;}
   if(b>f) {tmp=b; b=f; f=tmp;}
   if(c>g) {tmp=c; c=g; g=tmp;}
   if(d>h) {tmp=d; d=h; h=tmp;}
   if(a>c) {tmp=a; a=c; c=tmp;}
   if(b>d) {tmp=b; b=d; d=tmp;}
   if(a>b) {tmp=a; a=b; b=tmp;}
   if(e>g) {tmp=e; e=g; g=tmp;}
   if(f>h) {tmp=f; f=h; h=tmp;}
   if(g>h) {tmp=g; g=h; h=tmp;}
   if(a>i) {tmp=a; a=i; i=tmp;}
   if(h>i) {tmp=h; h=i; i=tmp;}
   return(med7(b, c, d, e, f, g, h));
}

/* Fast median */
template <typename DATA>
DATA median(int n, DATA *x)
{
   register int i, j, l;
   register DATA tmp;

/* Do insertion or heap median */
   if(n>7) goto HEAPMED;

/* Low n special cases (skip 8 and 9) */
   switch(n) {
      case 0:
         return(0.0);

      case 1:
         return(x[0]);

      case 2:
         return((x[0]+x[1])/2);

      case 3:
         return(med3(x[0], x[1], x[2]));

      case 4:
         return(med4(x[0], x[1], x[2], x[3]));

      case 5:
         return(med5(x[0], x[1], x[2], x[3], x[4]));

      case 6:
         return(med6(x[0], x[1], x[2], x[3], x[4], x[5]));

      case 7:
         return(med7(x[0], x[1], x[2], x[3], x[4], x[5], x[6]));

      case 8:
         return(med8(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7]));
         break;

      case 9:
         return(med9(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8]));
         break;

      default:
         ;;
         break;
   }

HEAPMED:
   if(n<INSHEAP) {

/* Insert-sort median */
      for(j=n-2; j>=0; j--) {
         i = n/2 + j + 1;
         if(i > n-1) i = n-1;
         for(l=j; l<i && x[j]>x[l+1]; l++);
         if(l > j) {
            tmp = x[j];
            for(i=j; i<l; i++) x[i] = x[i+1];
            x[l] = tmp;
         }
      }
      return( (n&1) ? x[n/2] : (x[n/2-1]+x[n/2])/2 );
   }

/* Default for large n is quantile() */
   return(quantile(n, x, NULL, 0.5));
}



////////////////////////////////////
// Typed versions of each routine //
////////////////////////////////////



/* quantile (e.g. median): 16 bit int data */
short int quantile_s(int n, short int *x, double *wgt, double q)
{
   return(quantile<short int>(n, x, wgt, q));
}

/* quantile (e.g. median): unsigned short data */
unsigned short quantile_u(int n, unsigned short *x, double *wgt, double q)
{
   return(quantile<unsigned short>(n, x, wgt, q));
}

/* quantile (e.g. median): int data */
int quantile_i(int n, int *x, double *wgt, double q)
{
   return(quantile<int>(n, x, wgt, q));
}

/* quantile (e.g. median): 64 bit int data */
long long int quantile_I(int n, long long int *x, double *wgt, double q)
{
   return(quantile<long long int>(n, x, wgt, q));
}

/* quantile (e.g. median): 32 bit float data */
float quantile_f(int n, float *x, double *wgt, double q)
{
   return(quantile<float>(n, x, wgt, q));
}

/* quantile (e.g. median): 64 bit double data */
double quantile_d(int n, double *x, double *wgt, double q)
{
   return(quantile<double>(n, x, wgt, q));
}

/* quantile (e.g. median): 128 bit long double data */
long double quantile_D(int n, long double *x, double *wgt, double q)
{
   return(quantile<long double>(n, x, wgt, q));
}




/* Quicksort, possibly carry an integer array: 16 bit int data */
void tsort_s(int n, short int *x, int *idx)
{
   qsort2<short int, int>(n, x, idx);
}

/* Quicksort, possibly carry a pointer array: 16 bit int data */
void tsort_sp(int n, short int *x, void* *idx)
{
   qsort2<short int, void *>(n, x, idx);
}

/* Quicksort, possibly carry an integer array: 16 bit unsigned int data */
void tsort_u(int n, unsigned short *x, int *idx)
{
   qsort2<unsigned short int, int>(n, x, idx);
}

/* Quicksort, possibly carry a pointer array: 16 bit unsigned int data */
void tsort_up(int n, unsigned short *x, void* *idx)
{
   qsort2<unsigned short int, void *>(n, x, idx);
}

/* Quicksort, possibly carry an integer array: integer data */
void tsort_i(int n, int *x, int *idx)
{
   qsort2<int, int>(n, x, idx);
}

/* Quicksort, possibly carry a pointer array: integer data */
void tsort_ip(int n, int *x, void* *idx)
{
   qsort2<int, void *>(n, x, idx);
}

/* Quicksort, possibly carry an integer array: 64 bit integer data */
void tsort_I(int n, long long int *x, int *idx)
{
   qsort2<long long int, int>(n, x, idx);
}

/* Quicksort, possibly carry a pointer array: 64 bit integer data */
void tsort_Ip(int n, long long int *x, void* *idx)
{
   qsort2<long long int, void *>(n, x, idx);
}

/* Quicksort, possibly carry an integer array: 32 bit float data */
void tsort_f(int n, float *x, int *idx)
{
   qsort2<float, int>(n, x, idx);
}

/* Quicksort, possibly carry a pointer array: 32 bit float data */
void tsort_fp(int n, float *x, void* *idx)
{
   qsort2<float, void *>(n, x, idx);
}

/* Quicksort, possibly carry an integer array: 64 bit double data */
void tsort_d(int n, double *x, int *idx)
{
   qsort2<double, int>(n, x, idx);
}

/* Quicksort, possibly carry a pointer array: 64 bit double data */
void tsort_dp(int n, double *x, void* *idx)
{
   qsort2<double, void *>(n, x, idx);
}

/* Quicksort, possibly carry an integer array: 128 bit double data */
void tsort_D(int n, long double *x, int *idx)
{
   qsort2<long double, int>(n, x, idx);
}

/* Quicksort, possibly carry a pointer array: 128 bit double data */
void tsort_Dp(int n, long double *x, void* *idx)
{
   qsort2<long double, void *>(n, x, idx);
}

/* Quicksort, possibly carry an integer array: character string data */
void tsort_str(int n, char **x, int *idx)
{
   qsort2_str<int>(n, x, idx);
}

/* Quicksort, possibly carry a pointer array: character string data */
void tsort_strp(int n, char **x, void* *idx)
{
   qsort2_str<void *>(n, x, idx);
}


/* Simple, fast median, short integer data */
short median_s(int n, short *x)
{
   return(median<short int>(n, x));
}

/* Simple, fast median, unsigned short integer data */
unsigned short median_u(int n, unsigned short *x)
{
   return(median<unsigned short int>(n, x));
}

/* Simple, fast median, integer data */
int median_i(int n, int *x)
{
   return(median<int>(n, x));
}

/* Simple, fast median, integer data */
int median_I(int n, long long int *x)
{
   return(median<long long int>(n, x));
}

/* Simple, fast median, float data */
float median_f(int n, float *x)
{
   return(median<float>(n, x));
}

/* Simple, fast median, double data */
double median_d(int n, double *x)
{
   return(median<double>(n, x));
}

/* Simple, fast median, long double data */
long double median_D(int n, long double *x)
{
   return(median<long double>(n, x));
}







// FORTRAN Wrappers
/* quantile (e.g. median): 16 bit int data */
short int quantile_s_(int *n, short int *x, double *q)
{
   return(quantile<short int>(*n, x, NULL, *q));
}

/* quantile (e.g. median): 16 bit int data */
short int wquantile_s_(int *n, short int *x, double *wgt, double *q)
{
   return(quantile<short int>(*n, x, wgt, *q));
}

/* quantile (e.g. median): int data */
int quantile_i_(int *n, int *x, double *q)
{
   return(quantile<int>(*n, x, NULL, *q));
}
/* quantile (e.g. median): int data */
int wquantile_i_(int *n, int *x, double *wgt, double *q)
{
   return(quantile<int>(*n, x, wgt, *q));
}

/* quantile (e.g. median): 32 bit float data */
float quantile_f_(int *n, float *x, double *q)
{
   return(quantile<float>(*n, x, NULL, *q));
}
/* quantile (e.g. median): 32 bit float data */
float wquantile_f_(int *n, float *x, double *wgt, double *q)
{
   return(quantile<float>(*n, x, wgt, *q));
}

/* quantile (e.g. median): 64 bit double data */
double quantile_d_(int *n, double *x, double *q)
{
   return(quantile<double>(*n, x, NULL, *q));
}
/* quantile (e.g. median): 64 bit double data */
double wquantile_d_(int *n, double *x, double *wgt, double *q)
{
   return(quantile<double>(*n, x, wgt, *q));
}

/* Quicksort, possibly carry an integer array: 16 bit int data */
void tsort_s_(int *n, short int *x)
{
   qsort2<short int, int>(*n, x, NULL);
}
void tsort_aux_s_(int *n, short int *x, int *idx)
{
   qsort2<short int, int>(*n, x, idx);
}

/* Quicksort, possibly carry an integer array: integer data */
void tsort_i_(int *n, int *x)
{
   qsort2<int, int>(*n, x, NULL);
}
/* Quicksort, possibly carry an integer array: integer data */
void tsort_aux_i_(int *n, int *x, int *idx)
{
   qsort2<int, int>(*n, x, idx);
}

/* Quicksort, possibly carry an integer array: 32 bit float data */
void tsort_f_(int *n, float *x)
{
   qsort2<float, int>(*n, x, NULL);
}
/* Quicksort, possibly carry an integer array: 32 bit float data */
void tsort_aux_f_(int *n, float *x, int *idx)
{
   qsort2<float, int>(*n, x, idx);
}

/* Quicksort, possibly carry an integer array: 64 bit double data */
void tsort_d_(int *n, double *x)
{
   qsort2<double, int>(*n, x, NULL);
}
void tsort_aux_d_(int *n, double *x, int *idx)
{
   qsort2<double, int>(*n, x, idx);
}


/* Simple, fast median, short integer data */
short median_s_(int n, short *x)
{
   return(median<short int>(n, x));
}
/* Simple, fast median, integer data */
int median_i_(int n, int *x)
{
   return(median<int>(n, x));
}
/* Simple, fast median, float data */
float median_f_(int n, float *x)
{
   return(median<float>(n, x));
}
/* Simple, fast median, double data */
double median_d_(int n, double *x)
{
   return(median<double>(n, x));
}




// Insertion sort wrappers, no real advantage over qsort2
#if 1
/* insertion sort, possibly carry an integer array, short data */
void insort2_s(int n, short *x, int *idx)
{
   insort2<short int, int>(n, x, idx);
}

/* insertion sort, possibly carry an integer array, ushort data */
void insort2_u(int n, unsigned short *x, int *idx)
{
   insort2<unsigned short int, int>(n, x, idx);
}

/* insertion sort, possibly carry an integer array, int data */
void insort2_i(int n, int *x, int *idx)
{
   insort2<int, int>(n, x, idx);
}

/* insertion sort, possibly carry an integer array, 64 bit int data */
void insort2_I(int n, long long int *x, int *idx)
{
   insort2<long long int, int>(n, x, idx);
}

/* insertion sort, possibly carry an integer array, 32 bit float data */
void insort2_f(int n, float *x, int *idx)
{
   insort2<float, int>(n, x, idx);
}

/* insertion sort, possibly carry an integer array, 64 bit double data */
void insort2_d(int n, double *x, int *idx)
{
   insort2<double, int>(n, x, idx);
}

/* insertion sort, possibly carry an integer array, 128 bit double data */
void insort2_D(int n, long double *x, int *idx)
{
   insort2<long double, int>(n, x, idx);
}
#endif


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


