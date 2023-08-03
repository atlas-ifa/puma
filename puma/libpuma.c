/* Fit a solar system orbit to data: Position Using Motion with Acceleration */
/* v1.21 - 200101 create a library */
/* v1.2  - 191212 modified to do a much better job on longer arcs */
/* v1.0  - 180605 John Tonry, simplified from ssorb.c */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "orbit.h"
#include "puma.h"

// #define UNWEIGHTED_XROCK	/* Define to use no weights in trajectory fit */

#define EARTH_MOON_TIDES   /* Define to treat Earth/Moon acc individually */

// #define USE_DIFSYS         /* provide real integrator functions? */

/* Global variables */
int PUMA_VERB=0;
int GNPT;		/* (global) number of data points */
DET *GDET;		/* (global) set of data points */
int GIDX;		/* (global) index of reference time */

// FIXME: Provide a way to allow user access to these?
int NDATE;		/* Fit estimated acc(t) with linear/quadratic 1/2 */
int PHITERM;		/* number of terms of phi(t) polynomial */
int ATERM;		/* number of terms of acceleration polynomial */
REAL VPRIOR=50;		/* [km/s] prior on velocity */
REAL RNGMIN=1e-4;	/* [AU] roll up chi^2 below this range */
REAL RNGMAX=3.0;	/* [AU] roll up chi^2 beyond this range */

int NXROCK=0;		// count of calls to xrock()
int DEBUG=0;

REAL LNRMIN=-10.0;	// Minimum ln(rAU)
REAL LNRMAX=+5.0;	// Maximum ln(rAU)
REAL LNRWID=7.5;	// Slope = (LNRMAX-LNRMIN)/2/LNRWID

/* Convert -inf < x < inf to lnr constrained to LNRMIN < ln(r) < LNRMAX */
#define LNRLIM(x) ( (LNRMAX-LNRMIN)/2 * \
		      (x)/SQRT(LNRWID*LNRWID+(x)*(x)) + \
		      (LNRMAX+LNRMIN)/2 )
/* Convert lnr constrained to LNRMIN < ln(r) < LNRMAX ln(r) to -inf < x < inf */
#define LIMLNR(x) ( LNRWID*(((x)-(LNRMAX+LNRMIN)/2)/(LNRMAX-LNRMIN)*2) / \
		    SQRT(1-(((x)-(LNRMAX+LNRMIN)/2)/(LNRMAX-LNRMIN)*2) * \
			   (((x)-(LNRMAX+LNRMIN)/2)/(LNRMAX-LNRMIN)*2)) )

// #define NEWFIT		/* Turn on restricted fit range */

#if 0
mongo

terminal
set \5 -10
set \6 5
set \4 4

def qq
set y x * x
set \0 \4 * \4
set y y + \0
set y y sqrt
set y x / y
set \0 \6 - \5
set \0 \0 / 2
set y y * \0
set \0 \6 + \5
set \0 \0 / 2
set y y + \0
end
 
def pp
set \0 \6 + \5
set \0 \0 /  2
set \1 \6 - \5
set \1 \1  /  2
set a y - \0
set a a  /  \1
set b a * a
set b 1 - b
set b b sqrt
set b a  /  b
set b b * \4
end
#endif

void xyzsun(REAL mjd, VEC *sun);
void planet(int k, double mjd, double *par, double *r, double *v,
	    double *hmag, double *albedo, char *name);
REAL obliquity(REAL mjd);

/* Rotate equatorial coords to ecliptic */
void eq2ecl(REAL mjd, REAL ra, REAL dec, REAL *lambda, REAL *beta);

/* Rotate ecliptic coords to equatorial */
void ecl2eq(REAL mjd, REAL lambda, REAL beta, REAL *ra, REAL *dec);

/* Rotate ecliptic coordinate vector to equatorial */
void ecl2eq_vec(REAL mjd, VEC *ECL, VEC *EQ);

/* Rotate ecliptic orbit state vector to equatorial */
void ecl2eq_orb(REAL mjd, ORBIT *ecl, ORBIT *eq);

/* Solar system / ecliptic positions for observer and observation */
void sslocation(DET *d);

/* Fit a great circle trajectory to N detections, return chi^2, g, v, and RMS */
REAL gctchi(DET *d0, int ndet, DET *d, REAL R, REAL V, VEC *g, REAL *xrms, REAL *trms);

/* Fit a great circle trajectory to N det, return angular velocity and pole */
void gctfit(DET *d0, int ndet, DET *d, REAL R, REAL V, VEC *g, REAL *phi);

/* Evaluate a great circle fit at a given MJD */
void gcteval(DET *d0, VEC *g, REAL *phi, DET *d, VEC *P);

/* Fit a trajectory to N detections, return fit params, RMS, and chi^2 */
REAL xrock(int ndet, DET *d, REAL R, REAL VR,
	   VEC *g, REAL *phi, VEC *X, REAL *xrms, REAL *trms);

/* Perturb the observed unit vector */
void xperturb(int ndet, DET *det, int doit);

/* Evaluate a given trajectory rooted on detection d0 on N detections */
void xeval(DET *d0, REAL rbary, REAL vrad, VEC *X, DET *d);

#ifdef USE_DIFSYS
void xintegrate(
   DET *d0,     /* detection at reference time */
   REAL r0,     /* [m] barycentric range at reference time */
   REAL vr0,    /* [m/s] barycentric radial velocity at reference time */
   VEC *X,      /* Position, velocity, and acceleration coefficients */
   DET *d);     /* detection time desired */
#endif

/* Definition of function that mini() minimizes */
typedef int MINIFUNC(int npar, REAL *par, int *usepar, REAL *value,
		     int doderiv, REAL *deriv, REAL *curve);

/* Quicksort long doubles */
// FIXME: REAL=long double?
void tsort_D(int n, REAL *x, int *idx);
/* Quicksort doubles */
void tsort_d(int n, double *x, int *idx);

/* Black box minimization using Levenberg-Marquard algorithm */
int mini(int npar, REAL *par, int *usepar,
         MINIFUNC func, REAL *cov, int nderiv, int *lambda, int verbose);

/* Gaussian random number */
REAL grand();

/* Diagonalize a 2x2 symmetric matrix */
/* (a b) = (c -s) (v1 0 ) ( c s)
 * (b d)   (s  c) (0  v2) (-s c)
 */
void diagsym(REAL a, REAL b, REAL d, REAL *v1, REAL *v2, REAL *phi);

/* Invert a 3x3 matrix in place */
int invert3(REAL *M);

/* Invert a 3x3 matrix in place */
int invert(int n, REAL *M, REAL *det);

/* Put a prior on VR of vsig ~ 0+/-50km/s */
void prioritize(REAL *cov, REAL *x0, REAL *v0, REAL vprior);

/* Calculate the chi^2 for a linearized estimate at given r, vr */
int xrockfunc(int npar, REAL *par, int *usepar, REAL *value,
	      int doderiv, REAL *deriv, REAL *curve);

/* Parse a CSV string with real variables, return array and number */
int parse_eph(char *arg, int maxvar, REAL *var);





/* Initialize puma, return a handle to the data and fit structure */
int puma_init(void **p)
{
   PUMA *puma;

/* Allocate a data and fit structure */
/* Note that calloc() zero's all coefficients */
   if( (puma = (PUMA *)calloc(1, sizeof(PUMA))) == NULL ) return(-1);

   puma->ndet = 0;
   puma->nalloc = PUMA_NALLOC;
   if( (puma->det = (DET *)calloc(PUMA_NALLOC, sizeof(DET))) == NULL ) return(-1);
   puma->ndate = 1;

   *p = (void *)puma;
   return(0);
}


/* Free the puma data and fit */
int puma_free(void *p)
{
   PUMA *puma=(PUMA *)p;

/* Allocate a data and fit structure */
   free(puma->det);
   free(puma);
   return(0);
}


/* Reset the data count within puma */
void puma_reset(void *p)
{
   PUMA *puma=(PUMA *)p;
   puma->ndet = puma->nfit = 0;
   return;
}


/* Provide an observation, return index and count */
int puma_obs(		// return 1-based index of this detection
   void *p,		// puma handle returned by puma_init()
   int eph,		// 0/1 for fitted data or ephemerides only
   double mjd,		// [day] MJD at instant of observation
   double ra,		// [deg] observed J2000 RA
   double dec,		// [deg] observed J2000 Dec
   double xerr,		// [arcsec] cross-track error estimate (0->eph req)
   double terr,		// [arcsec] tangential-track error estimate
   double lng,		// [deg] observatory (east) longitude
   double lat,		// [deg] observatory latitude
   double elev)		// [m] observatory elevation above geoid
{
   PUMA *puma=(PUMA *)p;
   int i, n=puma->ndet;
   REAL dr=atan(1)/45, MINGAP=0.4L, *buf;

/* Sanity checks */
   if(eph == 0 && xerr <= 0) {
      fprintf(stderr, "puma_obs called with xerr %10.3e <= 0\n", xerr);
      exit(1);
   }
   if(eph == 0 && terr <= 0) {
      fprintf(stderr, "puma_obs called with terr %10.3e <= 0\n", terr);
      exit(1);
   }


/* Ensure there is space */
   if(puma->ndet >= puma->nalloc) {
      puma->det = (DET *)realloc(puma->det, (PUMA_NALLOC+puma->nalloc)*sizeof(DET));
      if(puma->det == NULL) return(-1);
      puma->nalloc += PUMA_NALLOC;
   }

/* Add to the detection list */
   puma->det[n].obs.mjd = mjd;
   puma->det[n].obs.ra = ra * dr;	// [rad]
   puma->det[n].obs.dec = dec * dr;	// [rad]
   puma->det[n].obs.lng = lng * dr;	// [rad]
   puma->det[n].obs.lat = lat * dr;	// [rad]
   puma->det[n].obs.alt = elev;		// [m]
   puma->det[n].obs.psf = xerr;		// [arcsec] (overloaded usage)

   if(eph == 0) {
/* Uncertainty in cross-track position [rad-2] */
      puma->det[n].wgt = 1 / (xerr*xerr/3600/3600 * dr*dr);
/* obs.len => squared ratio of tangential/cross-track errors (overload) */
      puma->det[n].obs.len = terr*terr/xerr/xerr;
   } else {
      puma->det[n].obs.psf = puma->det[n].obs.len = puma->det[n].wgt = 0;
   }
   puma->det[n].obs.pa = 0;			// [deg]
   puma->det[n].obs.ut1utc = 0.0;		// [sec] (nominal)
   puma->det[n].obs.lst = lst(puma->det[n].obs.mjd, puma->det[n].obs.lng);

/* Precession from MJD to J2000 */
   precnutrot(puma->det[n].obs.mjd, puma->det[n].obs.ut1utc, puma->det[n].obs.pnr);

/* Fill in all the solar system information at this puma->det[n].obs.mjd */
   sslocation(puma->det+n);

   puma->ndet += 1;
   if(eph == 0) puma->nfit += 1;

// FIXME: allow user to specify reference detection
// FIXME: this assumes that the fit points all precede the eph points
/* Reference time */
   puma->id0 = puma->nfit / 2;


// FIXME: allow user to ndate
// FIXME: don't keep allocating and freeing!
/* Number of distinct dates */
   buf = (REAL *)calloc(puma->nfit, sizeof(REAL));
   for(i=0; i<puma->nfit; i++) {
      buf[i] = puma->det[i].obs.mjd - puma->det[puma->id0].obs.mjd;
   }
   tsort_D(puma->nfit, buf, NULL);

   puma->ndate = 1;
   for(i=1; i<puma->nfit; i++) {
      if(buf[i] - buf[i-1] > MINGAP) puma->ndate += 1;
   }
   free(buf);

// FIXME: globals, bleah...
   GIDX = puma->id0;
   GDET = puma->det;
   GNPT = puma->nfit;
   NDATE = puma->ndate;

   ATERM = 2;	// linear acceleration
   PHITERM = 2;	// linear phi(t)
   if(NDATE >= 3 && puma->nfit > 3) {  // 3+ dates, >3 obs
      ATERM = 3;	// quadratic acceleration
      PHITERM = 3;	// quadratic phi(t)
   }
   if(NDATE >= 4 && puma->nfit > 4) {  // 4+ dates, >4 obs
      PHITERM = 4;	// cubit phi(t)
   }
//   if(NDATE >= 5) {
//      PHITERM = 4;	// cubic phi(t)
//   }
   PHITERM = MIN(PUMA_MAXPHI, PHITERM);		// restrain overenthusiasm

   return(puma->ndet);
}


/* Save or retrieve the id field from the puma->det->obs structure */
void puma_id(void *p, int idx, int put, char *id)
{
   PUMA *puma=(PUMA *)p;

   if(put) {
      strncpy(puma->det[idx].obs.id, id, 64);
   } else {
      strcpy(id, puma->det[idx].obs.id);
   }
   return;
}



/* Dump the puma structure */
void puma_dump(
   void *p,		// puma handle returned by puma_init()
   FILE *fp)		// [stdout] output file descriptor
{
   PUMA *puma=(PUMA *)p;
   int i;
   REAL dr=ATAN(1)/45;

   if(fp == NULL) fp = stdout;

   fprintf(fp, "Data:\n");
   for(i=0; i<puma->ndet; i++) {
      fprintf(fp, "%3d %s %12.6Lf %9.5Lf %9.5Lf %8.4Lf %8.4Lf %8.0Lf %8.2Lf %8.3Lf\n", 
	     i, i<puma->nfit?"fit":"eph",
	     puma->det[i].obs.mjd,
	     puma->det[i].obs.ra/dr,
	     puma->det[i].obs.dec/dr,
	     puma->det[i].obs.lng/dr,
	     puma->det[i].obs.lat/dr,
	     puma->det[i].obs.alt,
	     puma->det[i].obs.psf,
	     puma->det[i].obs.len);
   }
			       
   fprintf(fp, "Control:  ");
   fprintf(fp, "ndet= %d nfit= %d id0= %d nalloc= %d ndate= %d\n",
	  puma->ndet, puma->nfit, puma->id0, puma->nalloc, puma->ndate);

   fprintf(fp, "Fit:  ");
   fprintf(fp, "t0= %12.6Lf R[AU]= %9.6Lf VR[km/s]= %8.3Lf\n",
	   puma->det[puma->id0].obs.mjd, puma->R/AU, puma->VR/1e3);

   for(i=0; i<5; i++) {
      fprintf(fp, "A%d=  %11.4Le %11.4Le %11.4Le\n",
	     i, puma->A[i].x, puma->A[i].y, puma->A[i].z);
   }
   for(i=0; i<2; i++) {
      fprintf(fp, "g%d=  %11.4Le %11.4Le %11.4Le\n",
	     i, puma->g[i].x, puma->g[i].y, puma->g[i].z);
   }
   fprintf(fp, "phi=");
   for(i=0; i<PHITERM; i++) fprintf(fp, " %11.4Le", puma->phi[i]);
   fprintf(fp, "\n");

   return;
}



/* Evaluate the best linear fit at a particular (ln(r),vr/r) */
double puma_chi(	// return chi^2/N
   void *p,		// puma handle returned by puma_init()
   double lnr,		// [ln(AU)] ln(r/AU) barycenter-object distance
   double wr,		// [/day] vr/r barycenter-object dln(r)/dt
   double rms[4])	// [","/day^2] crs/tng RMS [0,1], and curvatures [2,3]
{
   PUMA *puma=(PUMA *)p;
   double chi;
   REAL R, VR, xrms[2], trms[2];

   R = AU * exp(lnr);
   VR = wr * R / SECDAY;

/* Update fit parameters */
   puma->R = R;
   puma->VR = VR;

   chi = xrock(puma->nfit, puma->det, R, VR,
	       puma->g, puma->phi, puma->A, xrms, trms);

   rms[0] = xrms[0];
   rms[1] = trms[0];
   rms[2] = xrms[1];
   rms[3] = trms[1];

   return(chi);
}



/* Perturb fit and data according to their Gaussian uncertainties */
double puma_perturb(	// return chi^2/N
   void *p,		// puma handle returned by puma_init()
   double *lnr,		// [ln(AU)] ln(r/AU) barycenter-object distance
   double *wr,		// [/day] vr/r barycenter-object dln(r)/dt
   double dlnr,		// [ln(AU)] uncertainty in lnr
   double dwr,		// [/day] uncertainty in wr
   double rms[4])	// [","/day^2] crs/tng RMS [0,1], and curvatures [2,3]
{
   PUMA *puma=(PUMA *)p;
   double chi;
   REAL R, VR, xrms[2], trms[2];


/* Perturb the data */
   xperturb(puma->nfit, puma->det, 1);

/* Perturb the fit, this assumes negligible covariance! */
   *lnr += grand() * dlnr;
   *wr += grand() * dwr;

   R = AU * exp(*lnr);
   VR = (*wr) * R / SECDAY;

/* Get a new solution, use ndet to update solution for ephemeris points */
   chi = xrock(puma->ndet, puma->det, R, VR,
	       puma->g, puma->phi, puma->A, xrms, trms);

/* Restore the data */
   xperturb(puma->nfit, puma->det, 0);

   rms[0] = xrms[0];
   rms[1] = trms[0];
   rms[2] = xrms[1];
   rms[3] = trms[1];

   return(chi);
}


#ifdef NEWFIT

/* Execute a nlls fit over lnr,wr starting at (*lnr,*wr), return results */
double puma_fit(	// return chi^2/N
   void *p,		// puma handle returned by puma_init()
   double *lnr,		// [ln(AU)] ln(r/AU) barycenter-object distance
   double *wr,		// [/day] vr/r barycenter-object dln(r)/dt
   double *dlnr,	// [ln(AU)] ln(r/AU) uncertainty
   double *dwr,		// [/day] vr/r uncertainty
   double *rwcov)	// [-1:1] lnr-wr covariance
{
   PUMA *puma=(PUMA *)p;
   int i, usepar[NDIFEQ], lambda, miniverb, nderiv, err;
   double chi, dr=atan(1)/45;
   REAL R, VR, xrms[2], trms[2], par[NDIFEQ+1], cov[NDIFEQ*NDIFEQ];

/* Change to [m, m/s] */
   R = AU * exp(*lnr);			// [m]
   VR = (*wr) * dr / SECDAY * R;	// [m/s]

/* Initialize parameters */
   par[0] = *lnr;			// ln[AU]
   par[1] = VR/1e3 / (R/AU);		// [km/s/AU]

   usepar[0] = usepar[1] = 1;
   lambda = -2;
   miniverb = PUMA_VERB > 1;		/* Mini verbosity level */
   nderiv = 0;
/* do the fit */
   for(i=0; i<MAXMINITER; i++) {
      if( (err=mini(2, par, usepar, xrockfunc, cov, nderiv, 
		    &lambda, miniverb)) == 0 ) break;
   }
   printf("LAMBDA fit par,cov: %12.6Lf %12.6Lf %12.6Lf %12.6Lf %12.6Lf\n",
	  par[0], par[1], cov[0], cov[3], cov[2]);

/* Linear parameters */
   *lnr = LNRLIM(par[0]);		// [ln(AU)]
   R = AU * exp(*lnr);			// [m]
   VR = 1e3 * exp(*lnr) * par[1];	// [m/s]
   *wr = VR / R * SECDAY;		// [/day]

   printf("LAMBDA lnr wr R VR: %12.6f %12.6f %12.6Lf %12.6Lf\n",
	  *lnr, *wr, R, VR);

/* Uncertainties */
   *dlnr = cov[0];
   *dwr = cov[3];
   *rwcov = cov[1];

/* results: change units for wr error to [/day] */
   *dwr *= 1e3 / AU * SECDAY;

/* Update fit parameters */
   puma->R = R;				// [m]
   puma->VR = VR;			// [m/s]

   chi = xrock(puma->nfit, puma->det, puma->R, puma->VR,
	       puma->g, puma->phi, puma->A, xrms, trms);

   return(chi);
}

#else

/* Execute a nlls fit over lnr,wr starting at (*lnr,*wr), return results */
double puma_fit(	// return chi^2/N
   void *p,		// puma handle returned by puma_init()
   double *lnr,		// [ln(AU)] ln(r/AU) barycenter-object distance
   double *wr,		// [/day] vr/r barycenter-object dln(r)/dt
   double *dlnr,	// [ln(AU)] ln(r/AU) uncertainty
   double *dwr,		// [/day] vr/r uncertainty
   double *rwcov)	// [-1:1] lnr-wr covariance
{
   PUMA *puma=(PUMA *)p;
   int i, usepar[NDIFEQ], lambda, miniverb, nderiv, err;
   double chi, dr=atan(1)/45;
   REAL R, VR, xrms[2], trms[2], par[NDIFEQ+1], cov[NDIFEQ*NDIFEQ];

/* Change to [m, m/s] */
   R = AU * exp(*lnr);			// [m]
   VR = (*wr) * dr / SECDAY * R;	// [m/s]

/* Initialize parameters */
   par[0] = *lnr;			// ln[AU]
   par[1] = VR/1e3 / (R/AU);		// [km/s/AU]

   usepar[0] = usepar[1] = 1;
   lambda = -2;
   miniverb = PUMA_VERB > 1;		/* Mini verbosity level */
   nderiv = 0;
/* do the fit */
   for(i=0; i<MAXMINITER; i++) {
      if( (err=mini(2, par, usepar, xrockfunc, cov, nderiv, 
		    &lambda, miniverb)) == 0 ) break;
   }
//   printf("fitpar: %d %d %12.6Lf %12.6Lf %12.6Lf %12.6Lf %12.6Lf\n",
//	  err, lambda, par[0], par[1], cov[0], cov[3], cov[2]);

/* Sanity check */
   if(err > 0 ||
      par[0] > 1e6 || par[0] < -1e6 ||
      par[1] > 1e6 || par[1] < -1e6 ||
      cov[0] > 1e6 || cov[0] < -1e6 ||
      cov[1] > 1e6 || cov[1] < -1e6 ||
      cov[3] > 1e6 || cov[3] < -1e6) {

      par[0] = par[1] = cov[1] = 0;
      cov[0] = 100;
      cov[3] = 1e6;
   }


/* Linear parameters */
   R = AU * EXP(par[0]);			// [m]
   VR = 1e3 * EXP(par[0]) * par[1];		// [m/s]

/* results: change units */
   *lnr = par[0];			// [ln(AU)]
   *wr = VR / R * SECDAY;		// [/day]

/* Uncertainties */
   *dlnr = cov[0];
   *dwr = cov[3];
   *rwcov = cov[1];

/* results: change units for wr error to [/day] */
   *dwr *= 1e3 / AU * SECDAY;

/* Update fit parameters */
   puma->R = R;				// [m]
   puma->VR = VR;			// [m/s]

   chi = xrock(puma->nfit, puma->det, puma->R, puma->VR,
	       puma->g, puma->phi, puma->A, xrms, trms);

   return(chi);
}

#endif

/* Return information from detection array */
void puma_fitpar(
   void *p,		// puma handle returned by puma_init()
   int *nfit,		// number of fit points
   int *ndet,		// number of points, fit + ephemerides
   int *whichorigin,	// which detection is the origin (1:nfit)
   int *ndate,		// acceleration polynomial order
   double *t0,		// [MJD] time origin
   double X[3],		// [m] object fit location in SS, ecliptic coords
   double V[3],		// [m/s] object fit velocity in SS, ecliptic coords
   double A[9],		// [m/s^n] object fit acceleration in SS poly coeffs
   double gcpole[3],	// [] great circle pole unit vector
   double gcorigin[3],	// [] great circle origin unit vector
   double gcphi[PUMA_MAXPHI])	// [rad/s^n] time poly coeffs around GC
{
   PUMA *puma=(PUMA *)p;
   DET *d0=puma->det+puma->id0;

   *nfit = puma->nfit;
   *ndet = puma->ndet;
   *whichorigin = puma->id0 + 1;
   *ndate = puma->ndate;
   *t0 = d0->obs.mjd;

   X[0] = puma->A[0].x;
   X[1] = puma->A[0].y;
   X[2] = puma->A[0].z;

   V[0] = puma->A[1].x;
   V[1] = puma->A[1].y;
   V[2] = puma->A[1].z;

   A[0] = puma->A[1].x;
   A[1] = puma->A[1].y;
   A[2] = puma->A[1].z;
   A[3] = puma->A[2].x;
   A[4] = puma->A[2].y;
   A[5] = puma->A[2].z;
   A[6] = puma->A[3].x;
   A[7] = puma->A[3].y;
   A[8] = puma->A[3].z;

   gcpole[0] = puma->g[0].x;
   gcpole[1] = puma->g[0].y;
   gcpole[2] = puma->g[0].z;

   gcorigin[0] = puma->g[1].x;
   gcorigin[1] = puma->g[1].y;
   gcorigin[2] = puma->g[1].z;

   gcphi[0] = puma->phi[0];
   gcphi[1] = puma->phi[1];
   gcphi[2] = puma->phi[2];

   return;
}


/* Evaluate fit at a given time and return residuals */
double puma_resid(	// return sigma
   void *p,		// puma handle returned by puma_init()
   int which,		// evaluate which detection?
   double *rafit,	// [deg] RA from observer to object
   double *decfit,	// [deg] Dec from observer to object
   double *dist,	// [AU] distance from observer to object
   double *vrad,	// [km/s] radial velocity between observer and object
   double *wtan, 	// [deg/day] tangential velocity between obs and object
   double *vpa, 	// [deg] tang vel direction, E from N
   double *xerr,	// [arcsec] cross-track error
   double *terr,	// [arcsec] tangential-track error
   double *ragct,	// [deg] GCT RA from observer to object
   double *decgct)	// [deg] GCT Dec from observer to object
{
   PUMA *puma=(PUMA *)p;
   DET *det=puma->det+which-1, *d0=puma->det+puma->id0;
   double sigma;
   VEC P, U, V, W, THAT;
   REAL delta, rgct, vtot, vr, vra, vdec, rho, ra, dec, dr=ATAN(1)/45;

/* If an ephemeris request, need to compute the distance at det->obs.mjd */
   if(det->wgt == 0) xeval(d0, puma->R, puma->VR, puma->A, det);

/* Evaluate the great circle fit, return SS location */
   gcteval(d0, puma->g, puma->phi, det, &P);

/* GCT location wrt observatory */
   P.x -= det->O.x;
   P.y -= det->O.y;
   P.z -= det->O.z;
   rgct = SQRT(DOT(P, P));
   ecl2eq(det->obs.mjd, ATAN2(P.y, P.x), ASIN(P.z/rgct), &ra, &dec);

   *ragct = ra / dr;
   *decgct = dec / dr;

/* Fit location wrt observatory (ecliptic) */
   U.x = det->P.x - det->O.x;
   U.y = det->P.y - det->O.y;
   U.z = det->P.z - det->O.z;
   delta = SQRT(DOT(U, U));

/* RA,Dec of fit location from observer */
   ecl2eq(det->obs.mjd, ATAN2(U.y, U.x), ASIN(U.z/delta), &ra, &dec);

   *rafit = ra / dr;
   *decfit = dec / dr;
   *dist = delta / AU;
   
/* If a real observation, return fit residuals and sigma */
   if(det->wgt > 0) {
      *xerr = det->xresid*3600/dr;
      *terr = det->tresid*3600/dr;
      sigma = det->dsig;

/* Otherwise fake it with crs,tng GCT-fit differences and arcsec diff */
   } else {
      V.x = P.x/rgct - U.x/delta;	// difference of unit vectors: GCT-fit
      V.y = P.y/rgct - U.y/delta;
      V.z = P.z/rgct - U.z/delta;
      *xerr = DOT(V, puma->g[0]) / dr * 3600;
      *terr = TPROD(U, P, puma->g[0]) / delta/rgct / dr * 3600;
      sigma = DOT(P, U) / (rgct*delta);
      sigma = MAX(-1,MIN(1,sigma));
      sigma = ACOS(sigma) / dr * 3600;
   }

/* Velocities: object fit wrt observer (ecliptic) */
   V.x = det->V.x - det->U.x;
   V.y = det->V.y - det->U.y;
   V.z = det->V.z - det->U.z;
   vr = DOT(V,U) / delta;

   *vrad = vr / 1e3;
   *wtan = SQRT(DOT(V,V) - vr*vr) / delta / dr * SECDAY;

/* V is now ecliptic tangential velocity */
   V.x -= vr * U.x/delta;
   V.y -= vr * U.y/delta;
   V.z -= vr * U.z/delta;
   vtot = SQRT(DOT(V,V));

/* Unit vector to object in equatorial coordinates */
   W.x = U.x / delta;
   W.y = U.y / delta;
   W.z = U.z / delta;
   ecl2eq_vec(det->obs.mjd, &W, &W);

/* Theta_hat at this pointing */
   rho = SQRT(W.x*W.x+W.y*W.y);
   if(rho > 0) {
      THAT.x = -W.z * W.x/rho;
      THAT.y = -W.z * W.y/rho;
   } else {
      THAT.x = 1;
      THAT.y = 0;
   }
   THAT.z = rho;

/* Rotate ecliptic tangential velocity vector to equatorial */
   ecl2eq_vec(det->obs.mjd, &V, &W);

/* Dot it with theta_hat to get PA */
   *vpa = DOT(W,THAT)/vtot;
   *vpa = MAX(-1,MIN(1,*vpa));
   *vpa = ACOS(*vpa) / dr;

/* PA>0 if UxV has positive equatorial z. */
   P.x = XCRS(U, V);
   P.y = YCRS(U, V);
   P.z = ZCRS(U, V);
   vtot = SQRT(DOT(P,P));
   ecl2eq(det->obs.mjd, ATAN2(P.y, P.x), ASIN(P.z/vtot), &vra, &vdec);
   if(vdec > 0) *vpa =  ABS(*vpa);
   else         *vpa = -ABS(*vpa);

   return(sigma);
}


/* Return position and velocity with respect to the barycenter */
void puma_bary(
   void *p,		// puma handle returned by puma_init()
   int which,		// evaluate which detection?
   int integrate,       // [0/1] normal puma eval or true integration?
   double *mjd,	        // MJD of this obs
   double X[3],         // [m] cartesian object-barycenter displacement
   double V[3],         // [m/s] cartesian object-barycenter velocity
   double dX[3],	// [m] object-barycenter displacement
   double dV[3],	// [m/s] object-barycenter velocity
   double *sigma)       // [rad] point uncertainty
{
   int i;
   PUMA *puma=(PUMA *)p;
   DET *d=puma->det+which-1, *d0=puma->det+puma->id0;
   double sig=0;
//   double rho, vr, lat, phi, vth, vph, rxy, sig=0, pi=4*atan(1);
//   VEC U, TH, PH;

/* If an ephemeris request, need to compute the distance at det->obs.mjd */
   if(d->wgt == 0) {
#ifdef USE_DIFSYS
      if(integrate == 0) xeval(d0, puma->R, puma->VR, puma->A, d);
      else               xintegrate(d0, puma->R, puma->VR, puma->A, d);
#else
/* If an ephemeris request, need to compute the distance at det->obs.mjd */
      if(integrate == 1) fprintf(stderr, "puma integration not compiled\n");
      xeval(d0, puma->R, puma->VR, puma->A, d);
#endif
   }

   *mjd = d->obs.mjd;

/* Typical uncertainty in fit */
   for(i=0; i<puma->nfit; i++) sig += puma->det[i].wgt;
   *sigma = sqrt(puma->nfit/sig);

/* location wrt barycenter (ecliptic) */
   X[0] = d->P.x - d->E.x;
   X[1] = d->P.y - d->E.y;
   X[2] = d->P.z - d->E.z;

/* velocity wrt barycenter (ecliptic) */
   V[0] = d->V.x - d->W.x;
   V[1] = d->V.y - d->W.y;
   V[2] = d->V.z - d->W.z;

/* dX is observatory wrt barycenter */
   dX[0] = d->O.x - d->E.x;
   dX[1] = d->O.y - d->E.y;
   dX[2] = d->O.z - d->E.z;

/* dV is observatory wrt barycenter */
   dV[0] = d->U.x - d->W.x;
   dV[1] = d->U.y - d->W.y; 
   dV[2] = d->U.z - d->W.z; 

   return;
}


/* Return variance of unit vector wrt const for det with m0 <= mjd <= m1 */
double chiconst(int ndet, DET *d, double m0, double m1, int *n)
{
   int i;
   double wgt, w;
   VEC Xmean, Xvar;
   Xmean.x = Xmean.y = Xmean.z = Xvar.x = Xvar.y = Xvar.z = wgt = 0;
/* Compute chi^2/N for a constant location test */
   for(*n=i=0; i<ndet; i++) {
      if(d[i].obs.mjd < m0 || d[i].obs.mjd > m1) continue;
      w = d[i].wgt;
      wgt += w;
      Xmean.x += w * d[i].X.x;
      Xmean.y += w * d[i].X.y;
      Xmean.z += w * d[i].X.z;
      Xvar.x += w * d[i].X.x * d[i].X.x;
      Xvar.y += w * d[i].X.y * d[i].X.y;
      Xvar.z += w * d[i].X.z * d[i].X.z;
      *n += 1;
   }
   if(wgt <= 0) return(-1.0);
   if(*n < 2) return(0.0);
   Xvar.x -= Xmean.x*Xmean.x / wgt;
   Xvar.y -= Xmean.y*Xmean.y / wgt;
   Xvar.z -= Xmean.z*Xmean.z / wgt;
   return(Xvar.x + Xvar.y + Xvar.z);
}


/* Test whether points provided to puma are consistent with no motion */
double puma_const(	// return chi^2/N for all points at constant location
   void *p,		// puma handle returned by puma_init()
   double *chin)	// other chi^2/N variants: double const, ...
{
   PUMA *puma=(PUMA *)p;
   DET *d=puma->det;
   int i, n, n0, n1;
   double chin0, chin1;

/* Find the biggest time gap, mjdsep = middle of earlier and later group */
   double *buf = (double *)calloc(puma->nfit, sizeof(double));
   for(i=0; i<puma->nfit; i++) buf[i] = d[i].obs.mjd;
   tsort_d(puma->nfit, buf, NULL);

   double mjdsep=buf[0], gap=0;
   for(i=1; i<puma->nfit; i++) {
      if(buf[i] - buf[i-1] > gap) {
         mjdsep = 0.5 * (buf[i] + buf[i-1]);
         gap = buf[i] - buf[i-1];
      }
   }

/* First test all points at a common unit vector */
   chin0 = chiconst(puma->nfit, d, buf[0], buf[puma->nfit-1], &n);
   chin0 /= MAX(1, 2*n-2);

/* Next test for two constant positions across biggest time gap */
   chin1 = chiconst(puma->nfit, d, buf[0], mjdsep, &n0);
   chin1 += chiconst(puma->nfit, d, mjdsep, buf[puma->nfit-1], &n1);

   *chin = chin1 / MAX(1, 2*puma->nfit-4);

   free(buf);
   return(chin0);
}


/* Return selected Keplerian orbit parameters wrt Sun or Earth */
void puma_kepler(
   void *p,		// puma handle returned by puma_init()
   int earthctr,	// [0,1] wrt Sun (0) or Earth (1)
   double *E,		// [(km/s)^2] orbit energy
   double *ecc,		// [] orbit eccentricity
   double *incl,	// [deg] orbit inclination (ecliptic!)
   double *a,		// [AU,km] orbit semi-major axis
   double *period)	// [year,day] orbit period 
{   
   PUMA *puma=(PUMA *)p;
   DET *d0=puma->det+puma->id0;		// Evaluate at reference time
   ORBIT orb;
   REAL dr=ATAN(1)/45;

/* Earth barycenter orbit (ecliptic!) */
   if(earthctr) {
      orb.X.x = d0->P.x - d0->E.x;
      orb.X.y = d0->P.y - d0->E.y;
      orb.X.z = d0->P.z - d0->E.z;
      orb.V.x = d0->V.x - d0->W.x;
      orb.V.y = d0->V.y - d0->W.y;
      orb.V.z = d0->V.z - d0->W.z;
      orb.t0  = d0->obs.mjd - d0->obs.ltt/SECDAY;
      keparams(&orb, GMEARTH+GMMOON);
      *period = orb.period / SECDAY;	// [day] period
      *a = orb.a / 1e3;			// [km] semi-major axis

/* Solar system orbit */
   } else {
      orb.X.x = d0->P.x;
      orb.X.y = d0->P.y;
      orb.X.z = d0->P.z;
      orb.V.x = d0->V.x;
      orb.V.y = d0->V.y;
      orb.V.z = d0->V.z;
      orb.t0 = d0->obs.mjd - d0->obs.ltt/SECDAY;
      keparams(&orb, GMSUN);
      *period = orb.period / SECDAY / YEAR;	// [year] period
      *a = orb.a / AU;				// [AU] semi-major axis

   }

   *E = orb.E / 1e6;			// [(km/s)^2] energy
   *ecc = orb.ecc;			// [] eccentricity
   *incl = orb.incl / dr;		// [deg] inclination

   return;
}


/* Return information from detection array */
void puma_obsinfo(
   void *p,		// puma handle returned by puma_init()
   int which,		// evaluate which detection?
   int *eph,		// 0/1 for fitted data or ephemerides only
   double *t0,		// [day] origin MJD
   double *mjd,		// [day] MJD provided
   double *ra,		// [deg] RA provided or RA_GCT
   double *dec)		// [deg] Dec provided or Dec_GCT
{
   PUMA *puma=(PUMA *)p;
   DET *det=puma->det+which-1, *d0=puma->det+puma->id0;
   REAL dr=ATAN(1)/45;

   *eph = (det->obs.psf == 0);
   *t0 = d0->obs.mjd;
   *mjd = det->obs.mjd;
   *ra = det->obs.ra / dr;
   *dec = det->obs.dec / dr;
   return;
}


#define PERIGEE_LEN 100

/* Estimate parameters of closest approach */
double puma_perigee(	// [AU] return barycentric distance of closest approach
   void *p,		// puma handle returned by puma_init()
   double *mperi,	// [day] MJD of closest approach
   double *ra,		// [deg] RA of closest approach
   double *dec,		// [deg] Dec of closest approach
   double *vtan)	// [deg/day] angular velocity at closest approach
{
   PUMA *puma=(PUMA *)p;
   DET *det, *d0=(puma->det)+(puma->id0);
   int i, which, sign;
   double mjd, dist, dr=atan(1)/45;
   REAL tstep, *r, *dt, delta, vr, rafit, decfit;
   VEC U, V;

// Fit a non-zero parabola to the existing distances?


/* Allocate some space */
   r = (REAL *)calloc(PERIGEE_LEN, sizeof(REAL));
   dt = (REAL *)calloc(PERIGEE_LEN, sizeof(REAL));


   sign = (puma->VR > 0) ? +1 : -1;	// [-1,+1] for incoming, outgoing
   tstep = puma->R / 1e4 / SECDAY;	// [day] assumes incoming at 10km/s
   if(tstep > 2) tstep = 2;		// [day] max step = 2 day

/* Search for 2*PERIGEE_LEN days from t0 */
   for(i=0; i<PERIGEE_LEN; i++) {
      mjd = d0->obs.mjd - sign*i*tstep;

/* Create a test observation using the observatory location of the first obs */
      which = puma_obs(puma, 1, mjd, 0., 0., 0., 0.,
		    puma->det[0].obs.lng,
		    puma->det[0].obs.lat,
		    puma->det[0].obs.alt);

/* WARNING: puma-obs() might have moved puma->det via realloc()! */
      d0 = puma->det + puma->id0;

/* Get the location at det->obs.mjd and distance */
      det = puma->det + which - 1;
      xeval(d0, puma->R, puma->VR, puma->A, det);
      dt[i] = mjd - d0->obs.mjd;
      r[i] = det->rho;
//      printf("%3d %12.6f %8.3Lf %12.8Lf\n", i, mjd, dt[i], r[i]/AU);

/* Eliminate this test MJD */
      puma->ndet -= 1;
      if(i>1 && r[i] > r[i-1]) break;	// Ensure 3 points and outgoing
   }

/* Fit a parabola to the points that straddle the minimum */
   i = MIN(i-1, PERIGEE_LEN-2);
   mjd = sign*tstep/2*(r[i+1]-r[i-1])/(r[i+1]+r[i-1]-2*r[i]) + dt[i] + d0->obs.mjd;

   which = puma_obs(puma, 1, mjd, 0., 0., 0., 0.,
		    puma->det[0].obs.lng,
		    puma->det[0].obs.lat,
		    puma->det[0].obs.alt);

/* Get the location at det->obs.mjd and distance */
   det = puma->det + which - 1;
   xeval(d0, puma->R, puma->VR, puma->A, det);
   dt[i] = mjd - d0->obs.mjd;
   r[i] = det->rho;
//   printf("%3d %12.6f %8.3Lf %12.8Lf\n", i, mjd, dt[i], r[i]/AU);

/* Vector from observatory to object fit */
   U.x = det->P.x - det->O.x;
   U.y = det->P.y - det->O.y;
   U.z = det->P.z - det->O.z;
   delta = SQRT(DOT(U, U));

/* RA,Dec of fit location from observer */
   ecl2eq(det->obs.mjd, ATAN2(U.y, U.x), ASIN(U.z/delta), &rafit, &decfit);


   *mperi = mjd;
   dist = r[i]/AU;
   *ra = rafit / dr;
   *dec = decfit / dr;
   
/* Velocities: object wrt observer */
   V.x = det->V.x - det->U.x;
   V.y = det->V.y - det->U.y;
   V.z = det->V.z - det->U.z;
   vr = DOT(V,U) / delta;

   *vtan = SQRT(DOT(V,V) - vr*vr) / delta / dr * SECDAY;

/* Eliminate this test MJD */
   puma->ndet -= 1;

   free(r);
   free(dt);

   return(dist);
}







////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////



/* Given a barycenter dist and observer unit vector update bary unit vec */
REAL detbary(	// return = observer-det range
   REAL rbary,	// barycenter dist to detection
   DET *d)	// detection info
{
   REAL r, rdot;

/* Distance from observer to detection, solve quadratic */
   rdot = DOT(d->OB, d->X);
   r = rbary*rbary - d->rog*d->rog + rdot*rdot;
   if(r > 0) {
      r = SQRT(r) - rdot;
      if(r > 0) {
/* Unit vector from bary to rock */
	 d->XB.x = (d->OB.x + r * d->X.x) / rbary;
	 d->XB.y = (d->OB.y + r * d->X.y) / rbary;
	 d->XB.z = (d->OB.z + r * d->X.z) / rbary;

/* Return range to observer */
	 return(r);
      }
   }

/* Uh oh: rbary is inside of obs-barycenter distance -> garbage. */
   r = SQRT(DOT(d->OB, d->OB));
   d->XB.x = d->OB.x / r;
   d->XB.y = d->OB.y / r;
   d->XB.z = d->OB.z / r;
   return(r);
}



static REAL *T=NULL;

/* Fit a great circle trajectory to N det, return angular velocity and pole */
void gctfit(DET *d0,		// origin detection
	    int ndet,		// number of detections
	    DET *d,		// detection array
	    REAL R,		// [m] estimated distance
	    REAL V,		// [m/s] estimated radial velocity
	    VEC *g,		// 0=GC pole
	    REAL *phi)		// [rad/sec^n] angle [0], vel [1], accel [2]
{
   int i, j, k, n;
   REAL uxx, uyy, uzz, uxy, uxz, uyz, det;
   REAL r, dt, amin, rmin, tmin, rho, p, tpow, tpowk;
   REAL T0, T1, T2, T3, T4, q0, q1, q2;
   VEC w, U, b1;

// -O3 doesn't like such a large allocation for T[]
//   REAL T[PUMA_MAXPHI*PUMA_MAXPHI], Q[PUMA_MAXPHI];
   REAL Q[PUMA_MAXPHI];

   if(T == NULL) T = (REAL *)calloc(PUMA_MAXPHI*PUMA_MAXPHI, sizeof(REAL));

/* Quadratic distance: r(t) = amin*((t-t0)-tmin)^2 + RNGMIN */
   if(V != 0) {
//      rmin = RNGMIN*AU;
      rmin = 0.001 * R;
      tmin = -2*(R-rmin) / V;
      amin = V*V / (4*(R-rmin));
   } else {
      rmin = R;
      amin = tmin = 0;
   }

/* Great circle solution matrix */
   uxx = uxy = uxz = uyy = uyz = uzz = 0.0;


   for(i=n=0; i<ndet; i++) {
      if(d[i].wgt == 0) continue;

      dt = (d[i].obs.mjd - d0->obs.mjd) * SECDAY;

/* Linear estimate of rock distance wrt bary */
//      rho = R + V*dt;
//      printf("lin: %12.8Lf", rho/AU);

/* Non-negative quadratic estimate of barycentric range */
      rho = amin*(dt-tmin)*(dt-tmin) + rmin;
//      printf("  quad: %12.8Lf\n", rho/AU);


/* Observation location wrt barycenter */
      r = detbary(rho, d+i);

/* Unit vector from barycenter to object */
      U.x = d[i].XB.x;
      U.y = d[i].XB.y;
      U.z = d[i].XB.z;

/* Accumulate sums */
      uxx += U.x * U.x * d[i].wgt;
      uyy += U.y * U.y * d[i].wgt;
      uzz += U.z * U.z * d[i].wgt;
      uxy += U.x * U.y * d[i].wgt;
      uyz += U.y * U.z * d[i].wgt;
      uxz += U.x * U.z * d[i].wgt;

      n++;
   }

   if(n > 2) {
/* Pivot off of the smallest */
      if(ABS(uxx) < ABS(uyy) && ABS(uxx) < ABS(uzz)) {/* |uxx| is the smallest */
         det = uyy*uzz - uyz*uyz;
         w.x = 1.0;
         w.y = (uyz*uxz - uzz*uxy) / det;
         w.z = (uyz*uxy - uyy*uxz) / det;

      } else if(ABS(uyy) < ABS(uzz)) {                /* |uyy| is the smallest */
         det = uxx*uzz - uxz*uxz;
         w.x = (uxz*uyz - uzz*uxy) / det;
         w.y = 1.0;
         w.z = (uxz*uxy - uxx*uyz) / det;

      } else {                                        /* |uzz| is the smallest */
         det = uxx*uyy - uxy*uxy;
         w.x = (uxy*uyz - uyy*uxz) / det;
         w.y = (uxy*uxz - uxx*uyz) / det;
         w.z = 1.0;
      }      

   } else {
/* Simple cross product of first two */
      w.x = XCRS(d[0].XB, d[1].XB);
      w.y = YCRS(d[0].XB, d[1].XB);
      w.z = ZCRS(d[0].XB, d[1].XB);
   }

/* Normalize the great circle unit vector */
   det = SQRT(DOT(w,w));
   if(det != 0) {
      g[0].x = w.x / det;
      g[0].y = w.y / det;
      g[0].z = w.z / det;
   } else {
      g[0].x = g[0].y = 0;
      g[0].z = 1;
   }

/* Where are we?  Unit vector and observatory location */
   if(PUMA_VERB > 1) {
      printf("GC pole: %9.6Lf %9.6Lf %9.6Lf\n", g[0].x, g[0].y, g[0].z);
   }
   
#ifdef GCTFIT_NANDEBUG
   if( !(g[0].x >= 0 || g[0].x<0) ||
       !(g[0].y >= 0 || g[0].y<0) ||
       !(g[0].z >= 0 || g[0].z<0) ) {
      printf("Got a nan\n");
   }
   printf("  GC pole: %9.6Lf %9.6Lf %9.6Lf", g[0].x, g[0].y, g[0].z);
   printf("  det %10.3Le\n", det);
   printf("R %10.3Le  V %10.3Le", R, V);
   printf("  rmin %10.3Le tmin %10.3Le amin %10.3Le", rmin, tmin, amin);
   printf("  dt %10.3Le rho %10.3Le r %10.3Le\n", dt, rho, r);
   printf("u %10.3Le %10.3Le %10.3Le  %10.3Le %10.3Le %10.3Le\n", uxx, uyy, uzz, uxy, uyz, uzz);

   for(i=0; i<ndet; i++) {
      if(d[i].wgt == 0) {
         printf("%5d\n", i);
         continue;
      }
      dt = (d[i].obs.mjd - d0->obs.mjd) * SECDAY;
      rho = amin*(dt-tmin)*(dt-tmin) + rmin;
      r = detbary(rho, d+i);
      U.x = d[i].XB.x;
      U.y = d[i].XB.y;
      U.z = d[i].XB.z;
      printf("%5d %10.3Le %10.3Le %10.3Le %10.3Le %10.3Le %10.3Le\n",
             i, U.x, U.y, U.z, dt, rho, r);
   }
//      double *BAD=NULL;
//      r = *BAD;
#endif



/* b1 = normal wrt pole-d0 plane */
   b1.x = XCRS(g[0], d0->XB);
   b1.y = YCRS(g[0], d0->XB);
   b1.z = ZCRS(g[0], d0->XB);
   r = SQRT(DOT(b1, b1));
   b1.x /= r;
   b1.y /= r;
   b1.z /= r;

/* g[1] = projection of d0 onto g equator = origin at t=0 */
   g[1].x = XCRS(b1, g[0]);
   g[1].y = YCRS(b1, g[0]);
   g[1].z = ZCRS(b1, g[0]);

/* Linear fit */
   if(PHITERM <= 2) {
      T0 = T1 = T2 = q0 = q1 = 0.0;
   } else if(PHITERM == 3) {
      T0 = T1 = T2 = T3 = T4 = q0 = q1 = q2 = 0.0;
   } else {
      for(j=0; j<PHITERM; j++) {
	 Q[j] = 0;
	 for(k=0; k<PHITERM; k++) T[j+k*PHITERM] = 0;
      }
   }



   for(i=0; i<ndet; i++) {
      if(d[i].wgt == 0) continue;
/* Time of detection wrt t0 [sec] */
      dt = (d[i].obs.mjd - d0->obs.mjd) * SECDAY;

/* Normal wrt pole-d[i] plane */
      U.x = XCRS(g[0], d[i].XB);
      U.y = YCRS(g[0], d[i].XB);
      U.z = ZCRS(g[0], d[i].XB);
      r = SQRT(DOT(U, U));
      U.x /= r;
      U.y /= r;
      U.z /= r;

/* CCW angle around great circle wrt d0 */
      p = ATAN2(TPROD(b1, U, g[0]), DOT(b1, U));

/* Assemble basis function terms and data vector */
      if(PHITERM <= 2) {
	 T0 += d[i].wgt;
	 T1 += dt*d[i].wgt;
	 T2 += dt*dt*d[i].wgt;
	 q0 += p*d[i].wgt;
	 q1 += p*dt*d[i].wgt;
      } else if(PHITERM == 3) {
	 T0 += d[i].wgt;
	 T1 += dt*d[i].wgt;
	 T2 += dt*dt*d[i].wgt;
	 T3 += dt*dt*dt*d[i].wgt;
	 T4 += dt*dt*dt*dt*d[i].wgt;
	 q0 += p*d[i].wgt;
	 q1 += p*dt*d[i].wgt;
	 q2 += p*dt*dt*d[i].wgt;
      } else {
	 for(j=0, tpow=1; j<PHITERM; j++, tpow*=dt) {
	    Q[j] += p * tpow * d[i].wgt;
	    for(k=0, tpowk=tpow; k<=j; k++, tpowk*=dt) {
	       T[k+j*PHITERM] += tpowk * d[i].wgt;
	    }
	 }
      }
//      printf("%10.6Lf %10.6Lf\n", dt/SECDAY, p*180/3.14159);
   }

/* Solve for the polynomial coefficients */
   if(PHITERM <= 2) {
/* Expect wgt ~ (1")^-2 ~ 1e-11 */
      det = T0*T2 - T1*T1;
      phi[0] = (T2*q0 - T1*q1) / det;
      phi[1] = (T0*q1 - T1*q0) / det;
      phi[2] = 0;
   } else if(PHITERM == 3) {
      det = 2*T1*T2*T3 + T0*T2*T4 - T0*T3*T3 - T1*T1*T4 - T2*T2*T2;
      phi[0] = ((T2*T4-T3*T3)*q0 + (T2*T3-T1*T4)*q1 + (T1*T3-T2*T2)*q2)/det;
      phi[1] = ((T2*T3-T1*T4)*q0 + (T0*T4-T2*T2)*q1 + (T1*T2-T0*T3)*q2)/det;
      phi[2] = ((T1*T3-T2*T2)*q0 + (T1*T2-T0*T3)*q1 + (T0*T2-T1*T1)*q2)/det;
   } else {
/* Fill in basis function matrix and solve */
      for(j=0; j<PHITERM-1; j++)
	 for(k=j+1; k<PHITERM; k++) T[k+j*PHITERM] = T[j+k*PHITERM];
      i = invert(PHITERM, T, &det);
      for(j=0; j<PHITERM; j++)
	 for(k=0, phi[j]=0; k<PHITERM; k++) phi[j] += T[k+j*PHITERM] * Q[k];
   }

/* By convention the angular velocity according is RHR wrt the pole */
   if(phi[1] < 0) {
      for(j=0; j<PHITERM; j++) phi[j] *= -1;
      g[0].x *= -1;
      g[0].y *= -1;
      g[0].z *= -1;
   }   

/* Evaluate and save angular velocities for all detections */
   for(i=0; i<ndet; i++) {
/* Time of detection wrt t0 [sec] */
      dt = (d[i].obs.mjd - d0->obs.mjd) * SECDAY;
/* angular velocity */
      for(j=PHITERM-1, d[i].omega=0; j>=1; j--)
	 d[i].omega = dt*d[i].omega + j*phi[j];
//      for(j=PHITERM-1, p=0; j>=0; j--) p = dt*p + phi[j];
//	 printf("%10.6Lf %10.6Lf\n", dt/SECDAY, p*180/3.14159);
   }

   return;
}

/* Evaluate a great circle fit for a given detection */
void gcteval(DET *d0, VEC *g, REAL *phi, DET *d, VEC *P)
{
   int j;
   REAL p, S, C, dt;
   VEC w;

/* Time of detection wrt t0 [sec] */
   dt = (d->obs.mjd - d0->obs.mjd) * SECDAY;

/* Angle around the great circle */
//   p = phi[0] + dt*(phi[1]+dt*phi[2]);
   for(j=PHITERM-2, p=phi[PHITERM-1]; j>=0; j--) p = dt*p + phi[j];
   S = SIN(p);
   C = COS(p);

/* Orthogonal vector */
   w.x = XCRS(g[0], g[1]);
   w.y = YCRS(g[0], g[1]);
   w.z = ZCRS(g[0], g[1]);

/* Position in the solar system */
   P->x = d->E.x + d->rho * (S*w.x + C*g[1].x);
   P->y = d->E.y + d->rho * (S*w.y + C*g[1].y);
   P->z = d->E.z + d->rho * (S*w.z + C*g[1].z);

   return;
}




/* Fit a great circle trajectory to N detections, return chi^2, g, and RMS */
REAL gctchi(DET *d0,		// origin detection
	    int ndet,		// number of detections
	    DET *d,		// detection array
	    REAL R,		// [m] estimated distance
	    REAL V,		// [m/s] estimated radial velocity
	    VEC *g,		// 0=GC pole, 1=pt on GC at t0, 2=tang ang vel
	    REAL *xrms,		// [arcsec] cross RMS
	    REAL *trms)		// [arcsec] tangential RMS
{
   int i, j;
   REAL p, pobs, dt, chi, phi[PUMA_MAXPHI];
   REAL dr=ATAN(1.0)/45.0;

   gctfit(d0, ndet, d, R, V, g, phi);

/* chi is chi^2 for N points in two dimensions, 2*N-4 DOF */
   chi = 0.0;
/* Cross-track residuals from dot products between g and each vector */
   *xrms = *trms = 0.0;
   for(i=0; i<ndet; i++) {
/* Time of detection wrt t0 [sec] */
      dt = (d[i].obs.mjd - d0->obs.mjd) * SECDAY;

      d[i].xresid = ASIN(DOT(d[i].XB, g[0]));
      chi += d[i].xresid * d[i].xresid * d[i].wgt;
      *xrms += d[i].xresid * d[i].xresid;

      pobs = ATAN2(TPROD(d[i].XB, *g, d0->XB), DOT(d0->XB, d[i].XB));
      for(j=PUMA_MAXPHI-2, p=phi[PUMA_MAXPHI-1]; j>=0; j--) p = dt*p + phi[j];
//      d[i].tresid = p - (phi[0] + dt*(phi[1]+dt*phi[2]));
      d[i].tresid = pobs - p;
      chi += d[i].tresid * d[i].tresid * d[i].wgt * d[i].obs.len;
      *trms += d[i].tresid * d[i].tresid;
#ifdef JTEMP
      if(PUMA_VERB > -1) {
#else
      if(PUMA_VERB > 1) {
#endif
	 printf("%9.3Lf %9.3Lf  %12.1Lf %10.6Lf\n", d[i].xresid*3600/dr, d[i].tresid*3600/dr, dt, p/dr);
      }
   }

   if(PUMA_VERB > 1) {
//   if(PUMA_VERB > -1) {
      printf("%9.6Lf %12.4Le %12.4Le\n", phi[0]/dr, phi[1]/dr, phi[2]/dr);
   }
/* Evaluate cross-GC RMS [arcsec] */
   *xrms = SQRT(*xrms/ndet) * (3600/dr);
   *trms = SQRT(*trms/ndet) * (3600/dr);

   if(PUMA_VERB > 1) {
      printf("Chi: %9.3Lf  RMS: %9.3Lf %9.3Lf\n", chi, *xrms, *trms);
   }
   return(chi);
}



/* Given a det loc wrt bary return the differential grav acc wrt bary, det-bary */
REAL detgrav(	// return = grav component of range acceleration
   VEC *xbary,	// bary vector to detection
   DET *d,	// detection info
   VEC *acc)	// return grav acc wrt bary
{
   REAL rs, rg, rddot;
   VEC xsun;

/* det location wrt sun */
   xsun.x = xbary->x + d->E.x;
   xsun.y = xbary->y + d->E.y;
   xsun.z = xbary->z + d->E.z;

/* det acceleration wrt sun from sun and earth */
#ifdef EARTH_MOON_TIDES
   REAL rm;
   VEC xgeo, xmoon;

   xgeo.x = xsun.x - d->G.x;
   xgeo.y = xsun.y - d->G.y;
   xgeo.z = xsun.z - d->G.z;
   xmoon.x = xsun.x - d->M.x;
   xmoon.y = xsun.y - d->M.y;
   xmoon.z = xsun.z - d->M.z;

   rs = SQRT(DOT(xsun, xsun));
   rg = SQRT(DOT(xgeo, xgeo));
   rm = SQRT(DOT(xmoon, xmoon));

   acc->x = -GMSUN/(rs*rs*rs) * xsun.x - GMEARTH/(rg*rg*rg) * xgeo.x - GMMOON/(rm*rm*rm) * xmoon.x;
   acc->y = -GMSUN/(rs*rs*rs) * xsun.y - GMEARTH/(rg*rg*rg) * xgeo.y - GMMOON/(rm*rm*rm) * xmoon.y;
   acc->z = -GMSUN/(rs*rs*rs) * xsun.z - GMEARTH/(rg*rg*rg) * xgeo.z - GMMOON/(rm*rm*rm) * xmoon.z;

#else
   rs = SQRT(DOT(xsun, xsun));
   rg = SQRT(DOT(*xbary, *xbary));
   acc->x = -GMSUN/(rs*rs*rs) * xsun.x - GMEARTH/(rg*rg*rg) * xbary->x;
   acc->y = -GMSUN/(rs*rs*rs) * xsun.y - GMEARTH/(rg*rg*rg) * xbary->y;
   acc->z = -GMSUN/(rs*rs*rs) * xsun.z - GMEARTH/(rg*rg*rg) * xbary->z;
#endif



/* det acceleration wrt bary from sun and earth */
   acc->x -= d->A.x;
   acc->y -= d->A.y;
   acc->z -= d->A.z;

/* Radial component of acceleration wrt bary */
   rddot = DOT(*acc, *xbary) / rg;

   return(rddot);
}

/* Given a bary dist and vec for two obs, estimate a better second dist  */
REAL detdist(
   DET *d0,	// detection info for det0
   REAL r0,	// bary dist to det 0
   REAL r0dot,	// first derivative of bary dist to det 0
   DET *d1,	// detection info for det1
   REAL r1,	// bary dist to det 1
   REAL *r1dot)	// updated radial velocity at det1
{
   REAL rnew, gravacc, dr1dt2, dt;
   VEC xg1, a1;

/* Initial barycentric vector to det 1 */
   rnew = detbary(r1, d1);
   xg1.x = r1 * d1->XB.x;
   xg1.y = r1 * d1->XB.y;
   xg1.z = r1 * d1->XB.z;

   dt = (d1->obs.mjd-d0->obs.mjd) * SECDAY;

/* Gravitational acceleration wrt Earth at det 1 */
   gravacc = detgrav(&xg1, d1, &a1);
   dr1dt2 = gravacc + d1->omega*d1->omega*r1;

/* Save range acceleration */
   d1->rhoacc = dr1dt2;

// This is what it would cost to get a range accel derivative wrt range
// About 60usec per 2019XC iteration
#if 0
   REAL dd;
   VEC am, xgm;
   rnew = detbary(r1+100, d1);
   xgm.x = (r1+100) * d1->XB.x;
   xgm.y = (r1+100) * d1->XB.y;
   xgm.z = (r1+100) * d1->XB.z;

   rnew = detgrav(&xgm, d1, &am);
   dd = DOT(am, xgm)/(r1+100) + d1->omega*d1->omega*(r1+100);

   printf("%11.4Le\n", (dd-dr1dt2)/100);
#endif



//   printf("%12.6Lf %12.6Lf %12.6Lf\n", dt/SECDAY, d1->omega*d1->omega*r1*1e3, dr1dt2*1e3);


/* New estimate for r1 */
// Definition of r0dot as physical velocity instead of rate of change
// of time-delayed range.  Acceleration is so small that it doesn't
// change its effect significantly, acc term changes by:
// O(v*a*dt/c^2) ~ 30km/s*1e-5km/s^2*1e5s/(3e5km/s)^2 ~ 1e-10
   rnew = r0 + (r0dot * dt + (d0->rhoacc/3 + dr1dt2/6) * dt*dt) * (1+r0dot/cLIGHT);
//   rnew = r0 + (r0dot * dt + (r0ddot/3 + dr1dt2/6) * dt*dt);

/* New estimate of the physical radial velocity */
   *r1dot = r0dot + 2*(d0->rhoacc/3 + dr1dt2/6) * dt;

/* Return new estimate for r1 */
   return(rnew);
}




/* Get bary distance and obs range for each of ndet detections, linear acc */
void distances_acc1(DET *det0, int ndet, DET *d, REAL rbary, REAL vrad,
		    VEC *g, REAL *phi)
{
   int i;
   REAL dt, robs, gravacc, rho, vi, rmin, amin, tmin;
   VEC xg0, a0;

/* Location wrt bary x0 at reference det0, robs=observer range */
   robs = detbary(rbary, det0);
   xg0.x = rbary * det0->XB.x;
   xg0.y = rbary * det0->XB.y;
   xg0.z = rbary * det0->XB.z;

/* Gravitational acceleration at reference det */
   gravacc = detgrav(&xg0, det0, &a0);

   if(PUMA_VERB > 1) {
      printf("rbary,obs %12.8Lf %12.8Lf v[km/s] %8.4Lf\n",
	     rbary/AU, robs/AU, vrad/1e3);
   }

/* Acceleration of range wrt bary for det 0 and det 1 */
   det0->rhoacc = gravacc + det0->omega*det0->omega*rbary;

/* Quadratic distance: r(t) = amin*((t-t0)-tmin)^2 + RNGMIN */
   if(vrad != 0) {
//      rmin = RNGMIN*AU;
      rmin = 0.001*rbary;
      tmin = -2*(rbary-rmin) / vrad;
      amin = vrad*vrad / (4*(rbary-rmin));
   } else {
      rmin = rbary;
      tmin = amin = 0;
   }

//   printf("rbary,obs %12.8Lf %12.8Lf v[km/s] %8.4Lf\n",
//	  rbary/AU, robs/AU, vrad/1e3);

/* Estimate bary distance at each detection */
   for(i=0; i<ndet; i++) {

/* Time interval (sec) */
      dt = (d[i].obs.mjd-det0->obs.mjd) * SECDAY;

/* Linear estimate of rock distance wrt bary */
//      rho = rbary + vrad*dt;
//      printf(" lin:  %12.8Lf", rho/AU);

/* Non-negative quadratic estimate of barycentric range */
      rho = amin*(dt-tmin)*(dt-tmin) + rmin;
//      printf(" quad: %12.8Lf\n", rho/AU);
//      printf("%2d  %12.6Lf  %12.8Lf", i, dt/SECDAY, rho/AU);

/* Iterate once, applying an acceleration correction */
      rho = detdist(det0, rbary, vrad, d+i, rho, &vi);
//      printf("  %12.8Lf\n", rho/AU);

// This appears to matter at the ~1km level over a span of +/-2 days
/* Iterate once more, refining the acceleration correction */
#if 0
      printf("det %d:  %12.8Lf", i, d[i].obs.mjd, rho/AU);
      rho = detdist(det0, rbary, vrad, d+i, rho, &vi);
      printf("  %12.8Lf", rho/AU);
      rho = detdist(det0, rbary, vrad, d+i, rho, &vi);
      printf("  %12.8Lf\n", rho/AU);
#endif

/* Save the barycentric distance and observer range */
      d[i].rho = rho;
      d[i].r = detbary(rho, d+i);
   }
   
   return;
}

/* Get bary distance and obs range for each of ndet detections, quadratic acc */
void distances_acc2(DET *det0, int ndet, DET *d, REAL rbary, REAL vrad,
		    VEC *g, REAL *phi)
{
   int i;
   REAL T0, T1, T2, T3, T4, A0, A1, A2, Q0, Q1, Q2;
   REAL t0, ti, det, rho, acc;
   VEC U, S;

/* Fit a quadratic to the accelerations */
   t0 = det0->obs.mjd;
   T0 = T1 = T2 = T3 = T4 = Q0 = Q1 = Q2 = 0.0L;
   for(i=0; i<ndet; i++) {
      if(d[i].wgt <= 0) continue; 	/* skip eval only */
/* Add basis terms and data terms */
      ti = (d[i].obs.mjd-t0)*SECDAY - d[i].r/cLIGHT;

      Q0 += d[i].rhoacc;
      Q1 += ti * d[i].rhoacc;
      Q2 += ti*ti * d[i].rhoacc;

      T0 += 1;
      T1 += ti;
      T2 += ti*ti;
      T3 += ti*ti*ti;
      T4 += ti*ti*ti*ti;

//      printf("%10.6Lf %10.6Lf %10.6Lf\n", ti/SECDAY, d[i].rhoacc*1e3, d[i].rho/AU);
   }

/* Quadratic acceleration terms */
   det = 2*T1*T2*T3 + T0*T2*T4 - T0*T3*T3 - T1*T1*T4 - T2*T2*T2;
   A0 = ((T2*T4-T3*T3)*Q0 + (T2*T3-T1*T4)*Q1 + (T1*T3-T2*T2)*Q2)/det;
   A1 = ((T2*T3-T1*T4)*Q0 + (T0*T4-T2*T2)*Q1 + (T1*T2-T0*T3)*Q2)/det;
   A2 = ((T1*T3-T2*T2)*Q0 + (T1*T2-T0*T3)*Q1 + (T0*T2-T1*T1)*Q2)/det;
//   printf("%11.4Le %11.4Le %11.4Le\n", 1e3*A0, 1e3*A1*SECDAY, 1e3*A2*SECDAY*SECDAY);

/* Recompute ranges */
   for(i=0; i<ndet; i++) {
      if(d[i].wgt <= 0) continue; 	/* skip eval only */

      ti = (d[i].obs.mjd-t0)*SECDAY - d[i].r/cLIGHT;
      rho = rbary + ti*(vrad + ti/2*(A0 + ti/3*(A1 + ti/2*A2)));
      d[i].rho = rho;
      d[i].r = detbary(rho, d+i);

/* Update the accelerations */
      U.x = rho * d[i].XB.x;
      U.y = rho * d[i].XB.y;
      U.z = rho * d[i].XB.z;

/* Gravitational acceleration at reference det */
      acc = detgrav(&U, d+i, &S);

//      printf("%10.6Lf %10.6Lf %10.6Lf\n", 1e3*acc, 1e3*DOT(S,U)/rho, d[i].rho/AU);

/* Acceleration of range wrt bary */
      d[i].rhoacc = acc + d[i].omega*d[i].omega*rho;
   }

   return;
}


/* Get bary distance and obs range for each of ndet detections */
void distances(DET *det0, int ndet, DET *d, REAL rbary, REAL vrad,
	       VEC *g, REAL *phi)
{
/* GCT fit provides angular velocities wrt barycenter */
   gctfit(det0, ndet, d, rbary, vrad, g, phi);

/* Get bary distance and obs range for each of ndet detections, linear acc */
   if(NDATE <= 2) {
      distances_acc1(det0, ndet, d, rbary, vrad, g, phi);
   } else {
/* Get linear accelerations from acc1 */
      distances_acc1(det0, ndet, d, rbary, vrad, g, phi);
/* Iterate twice with a quadratic acceleration */
      distances_acc2(det0, ndet, d, rbary, vrad, g, phi);
      distances_acc2(det0, ndet, d, rbary, vrad, g, phi);
   }
   return;
}



/* Fit a trajectory to N detections, return chi^2/Ndof, X, V, and RMS */
REAL xrock(int ndet, DET *d, REAL rbary, REAL vrad,
	   VEC *g, REAL *phi, VEC *X, REAL *xrms, REAL *trms)
{
   int i, nuse;
   REAL t0, ti, dti, T0, T1, T2, T3, T4, r, chi, wgt, det, var;
   VEC X0, V0, A0, A1, A2;	// state vector and acc(ti)=A0+A1*ti+A2*ti^2
   VEC xfit, S, U, dU;		// temporary vectors
   REAL dr=ATAN(1.0)/45.0;

   NXROCK++;

/* Fill in bary distance and observer ranges for each detection */
   distances(d+GIDX, ndet, d, rbary, vrad, g, phi);


////////////////////////////////////////////////////////////////
/* Linear fit to acceleration */

/* time origin [day] */
   t0 = d[GIDX].obs.mjd;

   T0 = T1 = T2 = 0.0L;
   X0.x = X0.y = X0.z = V0.x = V0.y = V0.z = 0.0L;
   if(NDATE > 2) {
      T3 = T4 = U.x = U.y = U.z = 0.0L;
   } else {
      A2.x = A2.y = A2.z = 0;
   }
   wgt = 1.0;
   for(i=nuse=0; i<ndet; i++) {

      if(d[i].wgt <= 0) continue; 	/* skip eval only */

/* ti is sec wrt t0 but LTT in the past, so time of light emission */
      ti = (d[i].obs.mjd-t0)*SECDAY - d[i].r/cLIGHT;

// JTFIX
//      printf("Time: %11.3Le %11.3Le %11.3Le %11.3Le\n", d[i].obs.mjd, t0, d[i].r, ti);


/* Pretty accurate location in the solar system, given unit vector and range */
      xfit.x = d[i].O.x + d[i].r * d[i].X.x;
      xfit.y = d[i].O.y + d[i].r * d[i].X.y;
      xfit.z = d[i].O.z + d[i].r * d[i].X.z;
      r = SQRT(DOT(xfit, xfit));

/* Sun acceleration */
      A0.x = -GMSUN/r/r/r * xfit.x;
      A0.y = -GMSUN/r/r/r * xfit.y;
      A0.z = -GMSUN/r/r/r * xfit.z;

#ifdef EARTH_MOON_TIDES
/* Location of object with respect to Earth's center */
      xfit.x -= d[i].G.x;
      xfit.y -= d[i].G.y;
      xfit.z -= d[i].G.z;
      r = SQRT(DOT(xfit, xfit));

/* Sun and Earth acceleration */
      A0.x -= GMEARTH/r/r/r * xfit.x;
      A0.y -= GMEARTH/r/r/r * xfit.y;
      A0.z -= GMEARTH/r/r/r * xfit.z;

/* Location of object with respect to Moon's center */
      xfit.x -= d[i].M.x - d[i].G.x;
      xfit.y -= d[i].M.y - d[i].G.y;
      xfit.z -= d[i].M.z - d[i].G.z;
      r = SQRT(DOT(xfit, xfit));

/* Sun and Earth and Moon acceleration */
      A0.x -= GMMOON/r/r/r * xfit.x;
      A0.y -= GMMOON/r/r/r * xfit.y;
      A0.z -= GMMOON/r/r/r * xfit.z;
#else
/* Location of object with respect to Earth's center */
      xfit.x -= d[i].E.x;
      xfit.y -= d[i].E.y;
      xfit.z -= d[i].E.z;
      r = SQRT(DOT(xfit, xfit));

/* Sun and Earth acceleration */
      A0.x -= GMEARTH/r/r/r * xfit.x;
      A0.y -= GMEARTH/r/r/r * xfit.y;
      A0.z -= GMEARTH/r/r/r * xfit.z;
#endif

      if(PUMA_VERB == -1) {
	 printf("Data: %d %12.1Lf %12.8Lf %12.8Lf %12.8Lf\n", i, ti, A0.x*1e3, A0.y*1e3, A0.z*1e3);
      }

/* Add least squares components for all three axes */
      X0.x += A0.x * wgt;		// coopt X0,V0 temporarily
      X0.y += A0.y * wgt;
      X0.z += A0.z * wgt;
      V0.x += A0.x * ti * wgt;
      V0.y += A0.y * ti * wgt;
      V0.z += A0.z * ti * wgt;
      T0 += wgt;
      T1 += ti * wgt;
      T2 += ti*ti * wgt;
      if(NDATE > 2) {
	 U.x += A0.x * ti*ti * wgt;
	 U.y += A0.y * ti*ti * wgt;
	 U.z += A0.z * ti*ti * wgt;
	 T3 += ti*ti*ti * wgt;
	 T4 += ti*ti*ti*ti * wgt;
      }
// JTFIX
//      printf("%11.3Le %11.3Le %11.3Le\n", T0, T1, T2);
      nuse++;
   }

/* Mean linearized acceleration */

   if(NDATE > 2) {	/* Solve for quadratic acc(t) */
      det = 2*T1*T2*T3 + T0*T2*T4 - T0*T3*T3 - T1*T1*T4 - T2*T2*T2;
      A0.x = ((T2*T4-T3*T3)*X0.x + (T2*T3-T1*T4)*V0.x + (T1*T3-T2*T2)*U.x)/det;
      A1.x = ((T2*T3-T1*T4)*X0.x + (T0*T4-T2*T2)*V0.x + (T1*T2-T0*T3)*U.x)/det;
      A2.x = ((T1*T3-T2*T2)*X0.x + (T1*T2-T0*T3)*V0.x + (T0*T2-T1*T1)*U.x)/det;
      A0.y = ((T2*T4-T3*T3)*X0.y + (T2*T3-T1*T4)*V0.y + (T1*T3-T2*T2)*U.y)/det;
      A1.y = ((T2*T3-T1*T4)*X0.y + (T0*T4-T2*T2)*V0.y + (T1*T2-T0*T3)*U.y)/det;
      A2.y = ((T1*T3-T2*T2)*X0.y + (T1*T2-T0*T3)*V0.y + (T0*T2-T1*T1)*U.y)/det;
      A0.z = ((T2*T4-T3*T3)*X0.z + (T2*T3-T1*T4)*V0.z + (T1*T3-T2*T2)*U.z)/det;
      A1.z = ((T2*T3-T1*T4)*X0.z + (T0*T4-T2*T2)*V0.z + (T1*T2-T0*T3)*U.z)/det;
      A2.z = ((T1*T3-T2*T2)*X0.z + (T1*T2-T0*T3)*V0.z + (T0*T2-T1*T1)*U.z)/det;

      if(PUMA_VERB > 1 || PUMA_VERB == -1) {
	 printf("Accfit: %12.8Lf %12.8Lf %12.8Lf   %12.4Le %12.4Le %12.4Le   %12.4Le %12.4Le %12.4Le  %12.4Le\n",
		A0.x*1e3, A0.y*1e3, A0.z*1e3, 
		A1.x*1e3, A1.y*1e3, A1.z*1e3,
		A2.x*1e3, A2.y*1e3, A2.z*1e3, det);
      }

   } else {	/* Solve for linear acc(t) */
      det = T2*T0 - T1*T1;
      A0.x = (X0.x*T2 - V0.x*T1) / det;
      A1.x = (V0.x*T0 - X0.x*T1) / det;
      A0.y = (X0.y*T2 - V0.y*T1) / det;
      A1.y = (V0.y*T0 - X0.y*T1) / det;
      A0.z = (X0.z*T2 - V0.z*T1) / det;
      A1.z = (V0.z*T0 - X0.z*T1) / det;
      A2.x = A2.y = A2.z = 0;

// JTFIX
//      printf("%11.3Le %11.3Le %11.3Le\n", T0, T1, T2);

      if(PUMA_VERB > 1 || PUMA_VERB == -1) {
	 printf("Accfit: %12.8Lf %12.8Lf %12.8Lf   %12.4Le %12.4Le %12.4Le  %12.4Le\n",
		A0.x*1e3, A0.y*1e3, A0.z*1e3, 
		A1.x*1e3, A1.y*1e3, A1.z*1e3, det);
      }
   }
////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////
/* Fit linearized rock position with a state vector at t0-ltt */
   for(i=0; i<ndet; i++) {

      if(d[i].wgt <= 0) continue; 	/* skip eval only */

/* ti is sec wrt t0 but LTT in the past, so time of light emission */
      ti = (d[i].obs.mjd-t0)*SECDAY - d[i].r/cLIGHT;

/* Linearized location in the solar system */
      U.x = d[i].O.x + d[i].r*d[i].X.x - ti*ti/2*(A0.x+ti/3*(A1.x+ti*A2.x/2));
      U.y = d[i].O.y + d[i].r*d[i].X.y - ti*ti/2*(A0.y+ti/3*(A1.y+ti*A2.y/2));
      U.z = d[i].O.z + d[i].r*d[i].X.z - ti*ti/2*(A0.z+ti/3*(A1.z+ti*A2.z/2));

      if(PUMA_VERB > 1) {
	 printf("Rock: %9.1Lf  %9.6Lf %9.6Lf %9.6Lf\n",
		ti, U.x/AU, U.y/AU, U.z/AU);
      }

#ifndef UNWEIGHTED_X0	// default
      wgt = d[i].wgt;
#else
      wgt = 1.0;
#endif
      X0.x += U.x * wgt;
      X0.y += U.y * wgt;
      X0.z += U.z * wgt;
      V0.x += U.x * ti * wgt;
      V0.y += U.y * ti * wgt;
      V0.z += U.z * ti * wgt;
      T0 += wgt;
      T1 += ti * wgt;
      T2 += ti*ti * wgt;
   }
/* Mean linearized trajectory */
   det = T2*T0 - T1*T1;
   wgt = X0.x;
   X0.x = (X0.x*T2 - V0.x*T1) / det;
   V0.x = (V0.x*T0 - wgt*T1) / det;
   wgt = X0.y;
   X0.y = (X0.y*T2 - V0.y*T1) / det;
   V0.y = (V0.y*T0 - wgt*T1) / det;
   wgt = X0.z;
   X0.z = (X0.z*T2 - V0.z*T1) / det;
   V0.z = (V0.z*T0 - wgt*T1) / det;

   if(PUMA_VERB > 1) {
      printf("Rockloc: %9.6Lf %9.6Lf %9.6Lf   vel: %9.6Lf %9.6Lf %9.6Lf\n", 
	     X0.x/AU, X0.y/AU, X0.z/AU, 
	     V0.x*1e-3, V0.y*1e-3, V0.z*1e-3);
   }
////////////////////////////////////////////////////////////////


   

/*
  printf("Obs  loc: %9.6Lf %9.6Lf %9.6Lf   vel: %9.6Lf %9.6Lf %9.6Lf\n", 
  d[0].O.x/AU, d[0].O.y/AU, d[0].O.z/AU, 
  d[0].U.x*1e-3, d[0].U.y*1e-3, d[0].U.z*1e-3);
*/

////////////////////////////////////////////////////////////////
/* chi is chi^2 for N points in two dimensions, 2*N-4 DOF */
   chi = 0.0;

/* Cross-track and along-track residuals and curvature */
   xrms[0] = trms[0] = xrms[1] = trms[1] = T2 = 0.0;
   for(i=0; i<ndet; i++) {

/* ti is sec wrt t0 but LTT in the past, so time of light emission */
      d[i].obs.ltt = d[i].r/cLIGHT;
      dti = (d[i].obs.mjd-t0)*SECDAY;
      ti = dti - d[i].r/cLIGHT;

/* Save prediction of rock's location, including accelerations */
      d[i].P.x = X0.x + ti*(V0.x + ti/2*(A0.x + ti/3*(A1.x+ti*A2.x/2)));
      d[i].P.y = X0.y + ti*(V0.y + ti/2*(A0.y + ti/3*(A1.y+ti*A2.y/2)));
      d[i].P.z = X0.z + ti*(V0.z + ti/2*(A0.z + ti/3*(A1.z+ti*A2.z/2)));

/* Save the predicted velocity in the solar system */
      d[i].V.x = V0.x + ti*(A0.x + ti*(A1.x/2 + ti*A2.x/3));
      d[i].V.y = V0.y + ti*(A0.y + ti*(A1.y/2 + ti*A2.y/3));
      d[i].V.z = V0.z + ti*(A0.z + ti*(A1.z/2 + ti*A2.z/3));
      

//      printf("%.5Lf %.5Lf %.5Lf %.5Lf\n", ti, V0.x*ti/AU, V0.y*ti/AU, V0.z*ti/AU);
//      printf("%.5Lf %.5Lf %.5Lf\n", d[i].P.x/AU, d[i].P.y/AU, d[i].P.z/AU);

/* Predicted rock location wrt observer */
      S.x = d[i].P.x - d[i].O.x;
      S.y = d[i].P.y - d[i].O.y;
      S.z = d[i].P.z - d[i].O.z;
      r = SQRT(DOT(S, S));

/* Predicted unit vector from observatory to object */
      U.x = S.x / r;
      U.y = S.y / r;
      U.z = S.z / r;

//      printf("%.5Lf %.5Lf %.5Lf %.5Lf\n", ti, V0.x*ti/AU, V0.y*ti/AU, V0.z*ti/AU);
//      printf("%.5Lf %.5Lf %.5Lf\n", d[i].P.x/AU, d[i].P.y/AU, d[i].P.z/AU);

//      printf("%.5Lf %.5Lf %.5Lf\n", U.x, U.y, U.z);

/* Unit vector along relative velocity on sky */
      S.x = V0.x - d[i].U.x;		/* Velocity relative to observatory */
      S.y = V0.y - d[i].U.y;
      S.z = V0.z - d[i].U.z;
      r = DOT(S, U);
      S.x -= r * U.x;			/* Subtract radial component */
      S.y -= r * U.y;
      S.z -= r * U.z;
      r = SQRT(DOT(S, S));
      S.x /= r;				/* Renormalize */
      S.y /= r;
      S.z /= r;

/* Residual: observed unit vector minus prediction  */
      dU.x = d[i].X.x - U.x;
      dU.y = d[i].X.y - U.y;
      dU.z = d[i].X.z - U.z;

/* Cross-track residual, + is CCW wrt w */
      d[i].xresid = TPROD(U, dU, S);

/* Along-track residual, + is parallel to w */
      d[i].tresid = DOT(dU, S);

/* Location of object with respect to Earth's center */
      if(PUMA_VERB > 1) {
	 S.x = d[i].P.x - d[i].E.x;
	 S.y = d[i].P.y - d[i].E.y;
	 S.z = d[i].P.z - d[i].E.z;
	 r = SQRT(DOT(S, S));

	 printf("%12.3Lf %12.3Lf  %12.3Lf %12.3Lf\n",
		d[i].r/1e3, (d[i].r-r)/1e3, d[i].rho/1e3, (d[i].rho-r)/1e3);
      }
   
/* Accumulate contribution to chi^2 */
      wgt = 0;
      d[i].dsig = 0;	 	/* sigma error */
      if(d[i].wgt > 0) {
	 var = d[i].xresid * d[i].xresid * d[i].wgt;
	 chi += var;
	 d[i].dsig += var;
	 xrms[0] += d[i].xresid * d[i].xresid;

	 xrms[1] += d[i].xresid * dti*dti * d[i].wgt;
/* Recall d[i].obs.len has been coopted as tngerr*tngerr/(poserr*poserr) */
	 var = d[i].tresid * d[i].tresid * d[i].wgt / d[i].obs.len;
	 chi += var;
	 d[i].dsig += var;
	 trms[0] += d[i].tresid * d[i].tresid;
	 trms[1] += d[i].tresid * dti*dti * d[i].wgt / d[i].obs.len;

	 T2 += dti*dti*dti*dti * d[i].wgt;

	 wgt = DOT(dU, dU) * d[i].wgt;

	 d[i].use = 1;
	 d[i].dsig = SQRT(d[i].dsig);
      }
   }
////////////////////////////////////////////////////////////////


/* Evaluate cross-GC RMS [arcsec] */
   xrms[0] = SQRT(xrms[0]/nuse) * (3600/dr);
   trms[0] = SQRT(trms[0]/nuse) * (3600/dr);

/* Evaluate curvatures [arcsec/day^2] */
   xrms[1] *= 3600/dr / T2 * SECDAY*SECDAY;
   trms[1] *= 3600/dr / T2 * SECDAY*SECDAY;

/* Return state vector */
   X[0].x = X0.x;
   X[0].y = X0.y;
   X[0].z = X0.z;
   X[1].x = V0.x;
   X[1].y = V0.y;
   X[1].z = V0.z;
   X[2].x = A0.x;
   X[2].y = A0.y;
   X[2].z = A0.z;
   X[3].x = A1.x;
   X[3].y = A1.y;
   X[3].z = A1.z;
   X[4].x = A2.x;
   X[4].y = A2.y;
   X[4].z = A2.z;

/* Assume that this is part of a fit, eventually returning best r,vr */
   return(chi/MAX(1,(2*nuse-6)));
}

#ifdef NEWFIT

/* Calculate the chi^2 for a linearized estimate at given r, vr */
int xrockfunc(int npar, REAL *par, int *usepar, REAL *value,
	      int doderiv, REAL *deriv, REAL *curve)
{
   REAL xres[2], tres[2], phi[PUMA_MAXPHI], r;
   VEC X[5], g[2];
// par[0] = [AU] log(r)
// par[1] = [km/s/AU] vr/r

/* Stifle mini if it wants to blow up the exponential */
   r = exp(LNRLIM(par[0]));

   *value = xrock(GNPT, GDET, r*AU, par[1]*1e3*r, g, phi, X, xres, tres);

/* Apply a Gaussian prior of VPRIOR to the velocity likelihood */
   *value += par[1]*par[1]*r*r / (VPRIOR*VPRIOR);

//   printf("LAMBDA r,par[0],par[1]*r,limlnr: %11.6Lf %11.3Lf %11.3Lf %11.3Lf %12.4Le\n", r, par[0], par[1]*r, LIMLNR(log(r)), *value);

   return(0);
}

#else

/* Calculate the chi^2 for a linearized estimate at given r, vr */
int xrockfunc(int npar, REAL *par, int *usepar, REAL *value,
	      int doderiv, REAL *deriv, REAL *curve)
{
   REAL xres[2], tres[2], phi[PUMA_MAXPHI], r, vr;
   VEC X[5], g[2];
// par[0] = [AU] log(r)
// par[1] = [km/s/AU] vr/r

#if 0
/* Stifle mini if it wants to blow up the exponential */
   if(par[0] < 5) {
      r = exp(par[0]);
   } else {
      r = par[0]/5*exp(5.0);
   }
#endif

   r = exp(par[0]);
   if(r == 0) r = 1e-5;		// [AU]
   vr = 1e3*par[1] * r;
   *value = xrock(GNPT, GDET, r*AU, vr, g, phi, X, xres, tres);

//   printf("LAMBDA par[0],par[1]*r,r: %11.3Lf %11.3Lf %11.6Lf %12.4Le\n", par[0], par[1]*r, r, *value);

/* Apply a Gaussian prior of VPRIOR to the velocity likelihood */
   *value += par[1]*par[1]*r*r / (VPRIOR*VPRIOR);
   *value += par[1]*par[1]*r*r * par[1]*par[1]*r*r / (4*VPRIOR*VPRIOR*VPRIOR*VPRIOR);

/* Stifle r values beyond RNGMAX (3AU) */
   *value += r*r*r*r*r*r / (RNGMAX*RNGMAX*RNGMAX*RNGMAX*RNGMAX*RNGMAX);

/* Stifle r values inside of RNGMIN (1e-4AU = 2 Rearth)*/
   *value += (RNGMIN*RNGMIN*RNGMIN*RNGMIN*RNGMIN*RNGMIN) / (r*r*r*r*r*r);

//   printf("LAMBDA  r,par[1]*r: %11.3Lf %11.3Lf\n", r, par[1]*r);

   return(0);
}

#endif

/* Evaluate a given trajectory for a detection, origin detection d0 */
void xeval(DET *d0, REAL rbary, REAL vrad, VEC *X, DET *d)
{
   REAL ti;
   VEC U;

/* Initialize light travel times from that of d0 */
   d->obs.ltt = (rbary + (d->obs.mjd-d0->obs.mjd)*SECDAY*vrad) / cLIGHT;

/* Iterate once to get light travel time */
   ti = (d->obs.mjd-d0->obs.mjd)*SECDAY - d->obs.ltt;

/* Save prediction of rock's location, including accelerations */
   U.x = X[0].x + ti*(X[1].x + ti/2*(X[2].x + ti/3*(X[3].x+ti*X[4].x/2)));
   U.y = X[0].y + ti*(X[1].y + ti/2*(X[2].y + ti/3*(X[3].y+ti*X[4].y/2)));
   U.z = X[0].z + ti*(X[1].z + ti/2*(X[2].z + ti/3*(X[3].z+ti*X[4].z/2)));

/* ti is sec wrt t0 but LTT in the past, so time of light emission */
   U.x -= d->O.x;
   U.y -= d->O.y;
   U.z -= d->O.z;
   d->r = SQRT(DOT(U, U));
   ti = (d->obs.mjd-d0->obs.mjd)*SECDAY - d->r/cLIGHT;

/* Save prediction of rock's location, including accelerations */
   d->P.x = X[0].x + ti*(X[1].x + ti/2*(X[2].x + ti/3*(X[3].x+ti*X[4].x/2)));
   d->P.y = X[0].y + ti*(X[1].y + ti/2*(X[2].y + ti/3*(X[3].y+ti*X[4].y/2)));
   d->P.z = X[0].z + ti*(X[1].z + ti/2*(X[2].z + ti/3*(X[3].z+ti*X[4].z/2)));
   d->V.x = X[1].x + ti*(X[2].x + ti*(X[3].x/2 + ti*X[4].x/3));
   d->V.y = X[1].y + ti*(X[2].y + ti*(X[3].y/2 + ti*X[4].y/3));
   d->V.z = X[1].z + ti*(X[2].z + ti*(X[3].z/2 + ti*X[4].z/3));

/* Save barycentric distance, just in case */
   U.x = d->P.x - d->E.x;
   U.y = d->P.y - d->E.y;
   U.z = d->P.z - d->E.z;
   d->rho = SQRT(DOT(U, U));
   return;
}


/* Acceleration felt by object [m/s^2]   UTC is (*mjd+t/SECDAY) */
int accel(REAL t/*sec*/, REAL *P, REAL *Pdot, REAL *mjd/*start*/)
{
   double r, par[6], R[3], V[3], hmag, albedo;
   char name[80];
   VEC M, E;

/* Derivatives of first three variables */
   Pdot[0] = P[3];	/* dx/dt */
   Pdot[1] = P[4];	/* dy/dt */
   Pdot[2] = P[5];	/* dz/dt */

/* Acceleration due to Sun */
   r = SQRT(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
   Pdot[3] = -GMSUN * P[0] / (r*r*r);   /* d2x / dt2 */
   Pdot[4] = -GMSUN * P[1] / (r*r*r);   /* d2y / dt2 */
   Pdot[5] = -GMSUN * P[2] / (r*r*r);   /* d2z / dt2 */

/* Earth-Moon barycenter's position wrt Sun in ecliptic coords [AU] */
   planet(2, (double)(*mjd+t/SECDAY), par, R, V, &hmag, &albedo, name);

/* Moon position from center of Earth: ecliptic coords */
   lp_xyzmoon_ecliptic((double)(*mjd+t/SECDAY), &M);

/* Location of Earth wrt Sun in ecliptic coords */
   E.x = R[0]*AU - GMMOON/(GMEARTH+GMMOON) * M.x;
   E.y = R[1]*AU - GMMOON/(GMEARTH+GMMOON) * M.y;
   E.z = R[2]*AU - GMMOON/(GMEARTH+GMMOON) * M.z;

/* Location of Moon wrt Sun in ecliptic coords */
   M.x = R[0]*AU + GMEARTH/(GMEARTH+GMMOON) * M.x;
   M.y = R[1]*AU + GMEARTH/(GMEARTH+GMMOON) * M.y;
   M.z = R[2]*AU + GMEARTH/(GMEARTH+GMMOON) * M.z;

/* Acceleration due to Earth */
   r = SQRT((P[0]-E.x)*(P[0]-E.x)+(P[1]-E.y)*(P[1]-E.y)+(P[2]-E.z)*(P[2]-E.z));
   Pdot[3] -= GMEARTH * (P[0]-E.x) / (r*r*r);
   Pdot[4] -= GMEARTH * (P[1]-E.y) / (r*r*r);
   Pdot[5] -= GMEARTH * (P[2]-E.z) / (r*r*r);

/* Acceleration due to Moon */
   r = SQRT((P[0]-M.x)*(P[0]-M.x)+(P[1]-M.y)*(P[1]-M.y)+(P[2]-M.z)*(P[2]-M.z));
   Pdot[3] -= GMMOON * (P[0]-M.x) / (r*r*r);
   Pdot[4] -= GMMOON * (P[1]-M.y) / (r*r*r);
   Pdot[5] -= GMMOON * (P[2]-M.z) / (r*r*r);

   return(0);
}


#ifdef USE_DIFSYS

#define MAXITER 1000		/* Maximum number of difsys steps */
#define STEPMAX (50000.0)	/* [sec] Maximum initial difsys step */

/* Update X,V from initial conditions at mjd0 to mjdgoal */
int timestep(REAL mjd0, VEC *X, VEC *V, REAL mjdgoal)
{
   REAL w[NDIFEQ], ds, s, sgoal, eps=1e-12, h, ymax[NDIFEQ];
   int iter, err, error, newh;
//   printf(" %15.8Lf %15.8Lf\n", orb->t0, tgoal);
   if(mjdgoal == mjd0) return(0);
   w[0] = X->x;
   w[1] = X->y;
   w[2] = X->z;
   w[3] = V->x;
   w[4] = V->y;
   w[5] = V->z;
//   printf("%15.8Lf %10.1Lf %10.1Lf %10.1Lf   %10.3Lf %10.3Lf %10.3Lf\n", orb->t0, w[0], w[1], w[2], w[3], w[4], w[5]);
   s = 0.0;
   sgoal = (mjdgoal-mjd0) * SECDAY;
   h = STEPMAX;
   if(mjdgoal - mjd0 < 0) h = -STEPMAX;
   ymax[0] = ymax[1] = ymax[2] = ymax[3] = ymax[4] = ymax[5] = 0.0;
   for(iter=0; iter<MAXITER; iter++) {
      ds = sgoal - s;
      if(ABS(ds) > ABS(h)) ds = h;
//      printf("%10.1Lf %10.1Lf %10.5Lf %10.3Lf %10.3Lf %3d %3d", mjd0, sgoal, sgoal-s, h, ds, newh, iter);
      err = difsys(accel, NDIFEQ, &ds, &s, w, eps, ymax, &mjd0, &newh, &error);
//      printf("   %10.5Lf %10.1Lf %10.3Lf %10.3Lf %3d %3d\n", sgoal-s, s, h, ds, newh, iter);
      if(err && PUMA_VERB > 0) {
         fprintf(stderr, "difsys returned error %d at iter %d\n", err, iter);
      }
      h = ds;
      if(newh == 0 && ds > 0 && s >= sgoal) break;
      if(newh == 0 && ds < 0 && s <= sgoal) break;
   }
   if(iter == MAXITER) fprintf(stderr, "timestep: max iteration reached\n");
   if(err || ABS(sgoal-s) > 0.01) return(err);
   X->x = w[0];
   X->y = w[1];
   X->z = w[2];
   V->x = w[3];
   V->y = w[4];
   V->z = w[5];
//   printf("%10.1Lf %10.1Lf %10.1Lf   %10.3Lf %10.3Lf %10.3Lf\n", orb->X.x, orb->X.y, orb->X.z, orb->V.x, orb->V.y, orb->V.z);

   return(0);
}


/* Integrate a given trajectory for a detection, origin detection d0 */
void xintegrate(
   DET *d0,     /* detection at reference time */
   REAL r0,     /* [m] barycentric range at reference time */
   REAL vr0,    /* [m/s] barycentric radial velocity at reference time */
   VEC *X,      /* Position, velocity, and acceleration coefficients */
   DET *d)      /* detection time desired */
{
   int err;
   REAL t0, t1, ti, ltt;
   VEC P, V;

/* Light travel time from that of d0 */
   ltt = (r0 + (d->obs.mjd-d0->obs.mjd)*SECDAY*vr0) / cLIGHT;

/* Iterate once to get light travel time */
   ti = (d->obs.mjd-d0->obs.mjd)*SECDAY - ltt;

/* puma prediction of rock's location wrt Sun, including accelerations */
   P.x = X[0].x + ti*(X[1].x + ti/2*(X[2].x + ti/3*(X[3].x+ti*X[4].x/2)));
   P.y = X[0].y + ti*(X[1].y + ti/2*(X[2].y + ti/3*(X[3].y+ti*X[4].y/2)));
   P.z = X[0].z + ti*(X[1].z + ti/2*(X[2].z + ti/3*(X[3].z+ti*X[4].z/2)));

/* ti is sec wrt t0 but LTT in the past, so time of light emission */
   P.x -= d->O.x;
   P.y -= d->O.y;
   P.z -= d->O.z;
   d->r = SQRT(DOT(P, P));
   ti = (d->obs.mjd-d0->obs.mjd)*SECDAY - d->r/cLIGHT;

/* Now integrate properly from d0 to d */
   t0 = d0->obs.mjd - ltt/SECDAY;
   t1 = d->obs.mjd - d->r/cLIGHT/SECDAY;

/* Position wrt Sun at starting time t0; integrate to position at time t1 */
   P.x = d0->P.x;
   P.y = d0->P.y;
   P.z = d0->P.z;
   V.x = d0->V.x;
   V.y = d0->V.y;
   V.z = d0->V.z;

   if( (err=timestep(t0, &P, &V, t1)) ) {
      fprintf(stderr, "bad difsys %d\n", err);
   }

/* Save prediction of rock's location wrt Sun */
   d->P.x = P.x;
   d->P.y = P.y;
   d->P.z = P.z;
   d->V.x = V.x;
   d->V.y = V.y;
   d->V.z = V.z;

/* Save barycentric distance, just in case */
   P.x = d->P.x - d->E.x;
   P.y = d->P.y - d->E.y;
   P.z = d->P.z - d->E.z;
   d->rho = SQRT(DOT(P, P));
   return;
}
#endif

/* Gaussian random number */
REAL grand()
{
   int i=12;
   REAL g=-6.0*RAND_MAX;
   while(i--) g += random();
   return(g/RAND_MAX);
}


/* Diagonalize a 2x2 symmetric matrix */
/* (a b) = (c -s) (v1 0 ) ( c s)
 * (b d)   (s  c) (0  v2) (-s c)
 */
void diagsym(REAL a, REAL b, REAL d, REAL *v1, REAL *v2, REAL *phi)
{
   *phi = 0.5 * ATAN2(-2*b, a-d);
   *v1 = a*COS(*phi)*COS(*phi) - 2*b*SIN(*phi)*COS(*phi) + d*SIN(*phi)*SIN(*phi);
   *v2 = a*SIN(*phi)*SIN(*phi) + 2*b*SIN(*phi)*COS(*phi) + d*COS(*phi)*COS(*phi);
   return;
}

/* Put a prior on VR of vsig ~ 0+/-50km/s */
void prioritize(REAL *cov, REAL *x0, REAL *v0, REAL vsig)
{
   REAL cinv[4], new[4], det, a, b;

//   printf("cov\n%12.7Lf %12.7Lf\n%12.7Lf %12.7Lf\n", cov[0],cov[1],cov[2],cov[3]);

/* Inverse covariance matrix: assume cov in mini format of (sig1 r \ r sig2) */
   cinv[0] = 1/(1-cov[1]) / cov[0]/cov[0];
   cinv[1] = cinv[2] = -cov[1]/(1-cov[1]) / cov[0]/cov[3];
   cinv[3] = 1/(1-cov[1]) / cov[3]/cov[3];
//   printf("cinv\n%12.7Lf %12.7Lf\n%12.7Lf %12.7Lf\n", cinv[0],cinv[1],cinv[2],cinv[3]);

/* Modify with vsig, and invert to covariance matrix modified by v prior */
   det = cinv[0]*(cinv[3]+1/vsig/vsig) - cinv[1]*(cinv[2]);
   if(det == 0) {
      return;
   }
//   printf("det %12.7Lf\n", det);
   new[0] = (cinv[3]+1/vsig/vsig) / det;
   new[1] = new[2] = -cinv[1] / det;
   new[3] = cinv[0] / det;
//   printf("new\n%12.7Lf %12.7Lf\n%12.7Lf %12.7Lf\n", new[0],new[1],new[2],new[3]);

/* Compute new central values */
   a = cinv[0]*(*x0) + cinv[1]*(*v0);
   b = cinv[2]*(*x0) + cinv[3]*(*v0);

   *x0 = new[0]*a + new[1]*b;
   *v0 = new[2]*a + new[3]*b;

/* Update cov from new (in mini format)*/
   cov[0] = SQRT(ABS(new[0])) * (new[0]>0?+1:-1);
   cov[3] = SQRT(ABS(new[3])) * (new[0]>0?+1:-1);
   cov[1] = cov[2] = new[1] / cov[0] / cov[3];

   return;
}



/* Perturb the observed unit vector */
void xperturb(int ndet, DET *det, int doit)
{
   int i;
   REAL r, v, delta;
   VEC U, T, C;

   if(doit == 0) {
      for(i=0; i<ndet; i++) {
         if(det[i].xsave_valid) {       // safe to update X?
            det[i].X.x = det[i].XSAVE.x;
            det[i].X.y = det[i].XSAVE.y;
            det[i].X.z = det[i].XSAVE.z;
            det[i].xsave_valid = 0;     // X updated, XSAVE no longer safe
         }
      }
      return;
   }

   for(i=0; i<ndet; i++) {
      det[i].XSAVE.x = det[i].X.x;
      det[i].XSAVE.y = det[i].X.y;
      det[i].XSAVE.z = det[i].X.z;
      det[i].xsave_valid = 1;
      if(det[i].wgt <= 0) continue;

/* Fit location wrt observatory (ecliptic) */
      U.x = det->P.x - det->O.x;
      U.y = det->P.y - det->O.y;
      U.z = det->P.z - det->O.z;
      delta = SQRT(DOT(U, U));
      U.x /= delta;
      U.y /= delta;
      U.z /= delta;

/* Velocities: object fit wrt observer (ecliptic) */
      T.x = det->V.x - det->U.x;
      T.y = det->V.y - det->U.y;
      T.z = det->V.z - det->U.z;
      v = DOT(T,U);

/* T is now ecliptic tangential velocity */
      T.x -= v * U.x;
      T.y -= v * U.y;
      T.z -= v * U.z;

/* Divide to get tangential unit vector */
      v = SQRT(DOT(T,T));
      T.x /= v;
      T.y /= v;
      T.z /= v;
//   printf("%8.5Lf %8.5Lf %8.5Lf  %9.5f %9.5f  %8.5Lf %8.5Lf %8.5Lf  %9.5f %9.5f\n", U.x, U.y, U.z, atan2(det->X.y, det->X.x)*57.296, asin(det->X.z)*57.296, T.x, T.y, T.z, atan2(T.y, T.x)*57.296, asin(T.z)*57.296);

/* Cross product to get the cross unit vector */
      C.x = XCRS(U, T);
      C.y = YCRS(U, T);
      C.z = ZCRS(U, T);
//      printf("%8.5Lf %8.5Lf %8.5Lf  %8.5Lf %8.5Lf %8.5Lf  %8.5Lf %8.5Lf %8.5Lf\n", det[i].X.x, det[i].X.y, det[i].X.z, T.x, T.y, T.z, C.x, C.y, C.z);

/* Perturb the unit vector to the object from observer */
      det[i].X.x += C.x * SQRT(1/det[i].wgt) * grand();
      det[i].X.y += C.y * SQRT(1/det[i].wgt) * grand();
      det[i].X.z += C.z * SQRT(1/det[i].wgt) * grand();
      det[i].X.x += T.x * SQRT(det[i].obs.len/det[i].wgt) * grand();
      det[i].X.y += T.y * SQRT(det[i].obs.len/det[i].wgt) * grand();
      det[i].X.z += T.z * SQRT(det[i].obs.len/det[i].wgt) * grand();
      r = sqrt(DOT(det[i].X, det[i].X));
      det[i].X.x /= r;
      det[i].X.y /= r;
      det[i].X.z /= r;
   }
}

/* Rotate equatorial coords to ecliptic */
void eq2ecl(REAL mjd, REAL ra, REAL dec, REAL *lambda, REAL *beta)
{
   REAL C, S, pi=4*ATAN(1.0), incl;
   VEC X;

/* equatorial unit vector */
   X.x = COS(dec) * COS(ra);
   X.y = COS(dec) * SIN(ra);
   X.z = SIN(dec);

/* Rotate observation to ecliptic */
   incl = obliquity(mjd) * pi/180;
   C = X.y;
   S = X.z;
   X.y =  COS(incl) * C + SIN(incl) * S;
   X.z = -SIN(incl) * C + COS(incl) * S;

/* ecliptic coords */
   *lambda = ATAN2(X.y, X.x);
   if(*lambda < 0) *lambda += 2*pi;
   if(*lambda > 2*pi) *lambda -= 2*pi;
   *beta = ASIN(X.z);
   return;
}

/* Rotate ecliptic coords to equatorial */
void ecl2eq(REAL mjd, REAL lambda, REAL beta, REAL *ra, REAL *dec)
{
   REAL C, S, pi=4*ATAN(1.0), incl;
   VEC X;

/* equatorial unit vector */
   X.x = COS(beta) * COS(lambda);
   X.y = COS(beta) * SIN(lambda);
   X.z = SIN(beta);

/* Rotate observation to ecliptic */
   incl = obliquity(mjd) * pi/180;
   C = X.y;
   S = X.z;
   X.y =  COS(incl) * C - SIN(incl) * S;
   X.z =  SIN(incl) * C + COS(incl) * S;

/* ecliptic coords */
   *ra = ATAN2(X.y, X.x);
   if(*ra < 0) *ra += 2*pi;
   if(*ra > 2*pi) *ra -= 2*pi;
   *dec = ASIN(X.z);
   return;
}

/* Rotate ecliptic coordinate vector to equatorial */
void ecl2eq_vec(REAL mjd, VEC *ECL, VEC *EQ)
{
   REAL C, S, pi=4*ATAN(1.0), incl;

/* Rotate observation to ecliptic */
   incl = obliquity(mjd) * pi/180;
   C = (*ECL).y;
   S = (*ECL).z;
   (*EQ).x =  (*ECL).x;
   (*EQ).y =  COS(incl) * C - SIN(incl) * S;
   (*EQ).z =  SIN(incl) * C + COS(incl) * S;
   return;
}

/* Rotate ecliptic velocity to equatorial */
void vecl2eq(REAL mjd, REAL lambda, REAL beta, REAL *ra, REAL *dec)
{
   REAL C, S, pi=4*ATAN(1.0), incl;
   VEC X;

/* equatorial unit vector */
   X.x = COS(beta) * COS(lambda);
   X.y = COS(beta) * SIN(lambda);
   X.z = SIN(beta);

/* Rotate observation to ecliptic */
   incl = obliquity(mjd) * pi/180;
   C = X.y;
   S = X.z;
   X.y =  COS(incl) * C - SIN(incl) * S;
   X.z =  SIN(incl) * C + COS(incl) * S;

/* ecliptic coords */
   *ra = ATAN2(X.y, X.x);
   if(*ra < 0) *ra += 2*pi;
   if(*ra > 2*pi) *ra -= 2*pi;
   *dec = ASIN(X.z);
}

/* Rotate ecliptic orbit state vector to equatorial */
void ecl2eq_orb(REAL mjd, ORBIT *ecl, ORBIT *eq)
{
   REAL R, C, S, pi=4*ATAN(1.0), incl, beta, lambda;
   VEC X;

/* ecliptic unit vector */
   R = SQRT(ecl->X.x*ecl->X.x+ecl->X.y*ecl->X.y+ecl->X.z*ecl->X.z);
   X.x = ecl->X.x / R;
   X.y = ecl->X.y / R;
   X.z = ecl->X.z / R;
   beta = ASIN(X.z);
   lambda = ATAN2(X.y, X.x);

/* equatorial unit vector */
   X.x = COS(beta) * COS(lambda);
   X.y = COS(beta) * SIN(lambda);

/* Rotate observation to ecliptic */
   incl = obliquity(mjd) * pi/180;
   C = X.y;
   S = X.z;
   eq->X.x = R * X.x;
   eq->X.y = R * (COS(incl) * C - SIN(incl) * S);
   eq->X.z = R * (SIN(incl) * C + COS(incl) * S);

/* Rotate velocity to ecliptic */
   C = ecl->V.y;
   S = ecl->V.z;
   eq->V.x = ecl->V.x;
   eq->V.y = COS(incl) * C - SIN(incl) * S;
   eq->V.z = SIN(incl) * C + COS(incl) * S;

   eq->t0 = ecl->t0;
   return;
}

/* Solar system / ecliptic positions for observer and observation */
void sslocation(DET *d)
{
   REAL C, S, lam, x, y, z, omega, pi=4*ATAN(1.0), incl;
   VEC sun, moon;
   double par[6], r[3], v[3], hmag, albedo, dmjd;
   char name[80];

/* use unit vector for observation (equatorial) */
   d->X.x = COS(d->obs.ra) * COS(d->obs.dec);
   d->X.y = SIN(d->obs.ra) * COS(d->obs.dec);
   d->X.z = SIN(d->obs.dec);

/* Rotate observation to ecliptic */
   incl = obliquity(d->obs.mjd) * pi/180;
   C = d->X.y;
   S = d->X.z;
   d->X.y =  COS(incl) * C + SIN(incl) * S;
   d->X.z = -SIN(incl) * C + COS(incl) * S;

/* Change angular coords to ecliptic: "ra,dec" => lambda,beta */
   d->obs.lambda = ATAN2(d->X.y, d->X.x);
   if(d->obs.lambda < 0) d->obs.lambda += 2*pi;
   if(d->obs.lambda > 2*pi) d->obs.lambda -= 2*pi;
   d->obs.beta = ASIN(d->X.z);

/* compute observatory location in CIS*/

// From section K of the Astronomical Almanac:
// Geodetic coords:   lam, phi, h (height above geoid)
// Geocentric coords: lam, phi', rho
//   x = a rho cos(phi') cos(lam) = (aC+h) cos(phi) cos(lam)
//   y = a rho cos(phi') sin(lam) = (aC+h) cos(phi) sin(lam)
//   z = a rho sin(phi')          = (aS+h) sin(phi)
// where a is equatorial radius, C and S are functions of geodetic
// latitude phi and flattening f:
//   C = [cos^2(phi) + (1-f)^2 sin^2(phi)]^-1/2   S = (1-f)^2 C
// Polar radius b and eccentricity are given by:
//   b = a (1-f)  e^2 = 2f - f^2   1-e^2 = (1-f)^2
// For WGS84 a = 6378137m,  1/f=298.257223563,  GM=3.986005e14 m^3/s^2,
//          J2 = 0.00108263,  omega=7.292115e-5 rad/sec

   C = 1 / SQRT(COS(d->obs.lat)*COS(d->obs.lat) +
		(1-1/WGS84_finv)*(1-1/WGS84_finv) *
		SIN(d->obs.lat)*SIN(d->obs.lat));
   S = (1-1/WGS84_finv)*(1-1/WGS84_finv) * C;

/* Observatory location: Earth centered inertial */
   lam = d->obs.lng;
   x = (WGS84_a*C+d->obs.alt) * COS(lam) * COS(d->obs.lat);
   y = (WGS84_a*C+d->obs.alt) * SIN(lam) * COS(d->obs.lat);
   z = (WGS84_a*S+d->obs.alt) * SIN(d->obs.lat);

/* Geocentric observatory loc: precess and rotate to MJD in J2000 inertial */
   d->O.x = d->obs.pnr[0]*x + d->obs.pnr[1]*y + d->obs.pnr[2]*z;
   d->O.y = d->obs.pnr[3]*x + d->obs.pnr[4]*y + d->obs.pnr[5]*z;
   d->O.z = d->obs.pnr[6]*x + d->obs.pnr[7]*y + d->obs.pnr[8]*z;

/* Observatory rotational velocity, equatorial coords */
   omega = 2*pi / SECDAY * 1.0027379;
   d->U.x = -omega * d->O.y;
   d->U.y =  omega * d->O.x;
   d->U.z = 0.0;

/* Vector from Earth geocenter to Sun, equatorial coords */
//      xyzsun(d->obs.mjd, &sun);
   xyzsun_lpmoon(d->obs.mjd, &sun);

/* Vector from Sun to observatory, equatorial */
   d->O.x -= sun.x;
   d->O.y -= sun.y;
   d->O.z -= sun.z;

/* Vector from Sun to Earth geocenter, equatorial */
   d->E.x = -sun.x;
   d->E.y = -sun.y;
   d->E.z = -sun.z;

/* Rotate observatory location wrt sun to ecliptic */
   C = d->O.y;
   S = d->O.z;
   d->O.y =  COS(incl) * C + SIN(incl) * S;
   d->O.z = -SIN(incl) * C + COS(incl) * S;

/* Rotate geocenter location wrt sun to ecliptic */
   C = d->E.y;
   S = d->E.z;
   d->E.y =  COS(incl) * C + SIN(incl) * S;
   d->E.z = -SIN(incl) * C + COS(incl) * S;

/* Save geocenter wrt sun, ecl coords */
   d->G.x = d->E.x;
   d->G.y = d->E.y;
   d->G.z = d->E.z;

/* Vector from geo to observer */
   d->OB.x = d->O.x - d->E.x;
   d->OB.y = d->O.y - d->E.y;
   d->OB.z = d->O.z - d->E.z;
   d->rog = SQRT(DOT(d->OB, d->OB));

/* Earth-moon barycenter velocity wrt sun, ecliptic coords. v[AU/day] */
   dmjd = d->obs.mjd;
   planet(2, dmjd, par, r, v, &hmag, &albedo, name);
   v[0] *= AU/SECDAY;
   v[1] *= AU/SECDAY;
   v[2] *= AU/SECDAY;

   d->W.x = v[0];
   d->W.y = v[1];
   d->W.z = v[2];

/* Rotate observatory velocity wrt geocenter to ecliptic */
   C = d->U.y;
   S = d->U.z;
   d->U.y =  COS(incl) * C + SIN(incl) * S;
   d->U.z = -SIN(incl) * C + COS(incl) * S;

/* Add the barycenter ecliptic orbital velocity */
   d->U.x += v[0];
   d->U.y += v[1];
   d->U.z += v[2];

/* Earth-moon barycenter acceleration wrt sun, ecliptic coords. */
   dmjd = d->obs.mjd + 10.0/SECDAY;
   planet(2, dmjd, par, r, v, &hmag, &albedo, name);

   d->A.x = v[0]*AU/SECDAY;
   d->A.y = v[1]*AU/SECDAY;
   d->A.z = v[2]*AU/SECDAY;

   dmjd = d->obs.mjd - 10.0/SECDAY;
   planet(2, dmjd, par, r, v, &hmag, &albedo, name);

   d->A.x = (d->A.x - v[0]*AU/SECDAY) / 10 / 2;
   d->A.y = (d->A.y - v[1]*AU/SECDAY) / 10 / 2;
   d->A.z = (d->A.z - v[2]*AU/SECDAY) / 10 / 2;

/* Correct position and position-observatory to Earth-moon barycenter */
   lp_xyzmoon(d->obs.mjd, &moon);

#ifndef KEEP_GEOCENTER

/* Version that populates d->M and corrects(?) a bug */
#ifdef EARTH_MOON_TIDES
/* Say what?  Why am I not rotating from eq to ecl? */
   C = moon.y;
   S = moon.z;
   moon.y =  COS(incl) * C + SIN(incl) * S;
   moon.z = -SIN(incl) * C + COS(incl) * S;
/* Location of the moon wrt the Sun */
   d->M.x = d->G.x + moon.x;
   d->M.y = d->G.y + moon.y;
   d->M.z = d->G.z + moon.z;
#endif

   d->E.x += GMMOON/(GMMOON+GMEARTH) * moon.x;
   d->E.y += GMMOON/(GMMOON+GMEARTH) * moon.y;
   d->E.z += GMMOON/(GMMOON+GMEARTH) * moon.z;

   d->OB.x -= GMMOON/(GMMOON+GMEARTH) * moon.x;
   d->OB.y -= GMMOON/(GMMOON+GMEARTH) * moon.y;
   d->OB.z -= GMMOON/(GMMOON+GMEARTH) * moon.z;
   d->rog = SQRT(DOT(d->OB, d->OB));

#else
/* Correct velocity and acceleration of Earth-moon barycenter to geocenter */
/* Moon position wrt geocenter, equatorial coordinates */
//      xyzmoon(d->obs.mjd, &moon);
//      xyzmoon(d->obs.mjd+200.0/SECDAY, &moonv);
   lp_xyzmoon(d->obs.mjd+200.0/SECDAY, &moonv);
      
/* Rotate from equatorial to ecliptic */
   C = moon.y;
   S = moon.z;
   moon.y =  COS(incl) * C + SIN(incl) * S;
   moon.z = -SIN(incl) * C + COS(incl) * S;
   C = moonv.y;
   S = moonv.z;
   moonv.y =  COS(incl) * C + SIN(incl) * S;
   moonv.z = -SIN(incl) * C + COS(incl) * S;

/* Moon velocity wrt geocenter */
   moonv.x = (moonv.x - moon.x) / 200;
   moonv.y = (moonv.y - moon.y) / 200;
   moonv.z = (moonv.z - moon.z) / 200;

/* Include geo-barycenter velocity component in observatory velocity */
   d->U.x -= GMMOON/(GMMOON+GMEARTH) * moonv.x;
   d->U.y -= GMMOON/(GMMOON+GMEARTH) * moonv.y;
   d->U.z -= GMMOON/(GMMOON+GMEARTH) * moonv.z;

/* Include geo-barycenter velocity component in geocenter velocity */
   d->W.x -= GMMOON/(GMMOON+GMEARTH) * moonv.x;
   d->W.y -= GMMOON/(GMMOON+GMEARTH) * moonv.y;
   d->W.z -= GMMOON/(GMMOON+GMEARTH) * moonv.z;

/* Include geo-barycenter acceleration component in geocenter velocity */
   r2 = DOT(moon, moon);
   v2 = DOT(moonv, moonv);
   d->A.x += GMMOON/(GMMOON+GMEARTH) * v2/r2 * moon.x;
   d->A.y += GMMOON/(GMMOON+GMEARTH) * v2/r2 * moon.y;
   d->A.z += GMMOON/(GMMOON+GMEARTH) * v2/r2 * moon.z;
#endif
}

/* Solar system / ecliptic positions for observer and observation */
void cartesian(int npt, DET *d)
{
   int i;

   for(i=0; i<npt; i++) sslocation(d+i);

   if(PUMA_VERB > 1) {
      printf("    MJD         sec   xdet[]    y       z     Xgeo[AU]     y          z      Vx[km/s]      vy       vz   Ax[m/s2]     ay      az         a\n");
      for(i=0; i<npt; i++) {
	 printf("%11.5Lf  %6.0Lf  %7.4Lf %7.4Lf %7.4Lf  %9.5Lf %9.5Lf %9.5Lf  %8.3Lf %8.3Lf %8.3Lf  %8.5Lf %8.5Lf %8.5Lf %11.9Lf\n", 
		d[i].obs.mjd, SECDAY*(d[i].obs.mjd-d[0].obs.mjd),
		d[i].X.x, d[i].X.y, d[i].X.z, 
		d[i].E.x/AU, d[i].E.y/AU, d[i].E.z/AU,
		1e-3*d[i].U.x, 1e-3*d[i].U.y, 1e-3*d[i].U.z,
		d[i].A.x, d[i].A.y, d[i].A.z,
		SQRT(d[i].A.x*d[i].A.x+d[i].A.y*d[i].A.y) );
      }
   }
}

/* Invert a 3x3 matrix in place */
int invert3(REAL *M)
{
   REAL tprod;
   REAL A[3]={M[0],M[1],M[2]};
   REAL B[3]={M[3],M[4],M[5]};
   REAL C[3]={M[6],M[7],M[8]};

   tprod = (A[1]*B[2]-A[2]*B[1]) * C[0] +
           (A[2]*B[0]-A[0]*B[2]) * C[1] +
           (A[0]*B[1]-A[1]*B[0]) * C[2];

   if(tprod == 0) return(-1);

   M[0] = (B[1]*C[2]-B[2]*C[1]) / tprod;
   M[1] = (B[2]*C[0]-B[0]*C[2]) / tprod;
   M[2] = (B[0]*C[1]-B[1]*C[0]) / tprod;
   
   M[3] = (C[1]*A[2]-C[2]*A[1]) / tprod;
   M[4] = (C[2]*A[0]-C[0]*A[2]) / tprod;
   M[5] = (C[0]*A[1]-C[1]*A[0]) / tprod;
   
   M[6] = (A[1]*B[2]-A[2]*B[1]) / tprod;
   M[7] = (A[2]*B[0]-A[0]*B[2]) / tprod;
   M[8] = (A[0]*B[1]-A[1]*B[0]) / tprod;
   
   return(0);
}
