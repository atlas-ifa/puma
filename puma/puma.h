/* Header file for libpuma.c */


/* Usage of libpuma

The normal usage of libpuma is 

 1. Create a puma structure with puma_init()

 2. Provide libpuma with data via puma_obs().  Note that this is the
    avenue for both observational data which must be provided FIRST,
    and times for ephemeris requests.

 3. Evaluate goodness of fit at a particular (r,vr) with puma_chi()

 4. Execute a fit with puma_fit()

 5. Obtain desired results:

    Get fit parameters with puma_fitpar()

    Get observation residuals with puma_resid()

    Get ephemeris predictions with puma_resid()

    Get fit parameters with puma_fitpar()

    Explore statistical possibilities at an ephemeris MJD by perturbing
    the data and fit with puma_perturb()

    Get an estimate of closest approach with puma_perigee()

    Get Keplerian orbit parameters of the fit with puma_keper()

 6. Free the puma data with puma_free()

 */


// NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES NOTES

// The observer time origin is defined by det[id0].obs.mjd.  The fit assumes
// a barycentric range and range-dot of R,VR at that observer time.
//
//
// The fitted solar system position of the object is given by A[5], note
// that the position, velocity, acceleration, and jerk are all present,
// but that A[4] depends on ndate==3, and it is the quadratic term of
// a(t), without the 1/2 if it were the 2nd time derivative of a(t).
// t is light travel time lagged relative to the observed t'.

// The great circle fit has its pole g[0] and origin g[1].  The location
// around the great circle is given by cos(phi)*g[1]+sin(phi)*g[0]xg[1],
// where phi(t') = phi[0]+dt'*phi[1]+dt'*dt'*phi[2].  If ndate==2 phi[2]=0.





// FIXME: this must go into the puma structure in orbit.h, bad!
#define PUMA_MAXPHI 4	// Maximum number of phi terms for GCt


/* Initialize puma, return a handle to the data and fit structure */
int puma_init(void **p);

/* Free the puma data and fit */
int puma_free(void *p);

/* Reset the data count within puma */
void puma_reset(void *p);

/* Dump the puma structure to output file descriptor fp [NULL -> stdout] */
void puma_dump(void *p, FILE *fp);


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
   double elev);	// [m] observatory elevation above geoid


/* Evaluate the best linear fit at a particular (ln(r),vr/r) */
double puma_chi(	// return chi^2/N
   void *p,		// puma handle returned by puma_init()
   double lnr,		// [ln(AU)] ln(r/AU) barycenter-object distance
   double wr,		// [/day] vr/r barycenter-object dln(r)/dt
   double rms[4]);	// [","/day^2] crs/tng RMS [0,1], and curvatures [2,3]


/* Execute a nlls fit over lnr,wr starting at (*lnr,*wr), return results */
double puma_fit(	// return chi^2/N
   void *p,		// puma handle returned by puma_init()
   double *lnr,		// [ln(AU)] ln(r/AU) barycenter-object distance
   double *wr,		// [/day] vr/r barycenter-object dln(r)/dt
   double *dlnr,	// [ln(AU)] ln(r/AU) uncertainty
   double *dwr,		// [/day] vr/r uncertainty
   double *rwcov);	// [-1:1] lnr-wr covariance


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
   double gcphi[PUMA_MAXPHI]);	// [rad/s^n] time poly coeffs around GC


/* Perturb fit and data according to their Gaussian uncertainties */
double puma_perturb(	// return chi^2/N
   void *p,		// puma handle returned by puma_init()
   double *lnr,		// [ln(AU)] ln(r/AU) barycenter-object distance
   double *wr,		// [/day] vr/r barycenter-object dln(r)/dt
   double dlnr,		// [ln(AU)] uncertainty in lnr
   double dwr,		// [/day] uncertainty in wr
   double rms[4]);	// [","/day^2] crs/tng RMS [0,1], and curvatures [2,3]


/* Test whether points provided to puma are consistent with no motion */
double puma_const(	// return chi^2/N for all points at constant location
   void *p,		// puma handle returned by puma_init()
   double *chin);	// other chi^2/N variants: double const, ...


/* Evaluate fit for a particular obs/eph and return residuals */
double puma_resid(	// return sigma
   void *p,		// puma handle returned by puma_init()
   int which,		// detection (1:ndet) to evaluate
   double *rafit,	// [deg] RA from observer to object
   double *decfit,	// [deg] Dec from observer to object
   double *dist,	// [AU] distance from observer to object
   double *vrad,	// [km/s] radial velocity between observer and object
   double *wtan, 	// [deg/day] tangential velocity between obs and object
   double *vpa, 	// [deg] tang vel direction, E from N
   double *xerr,	// [arcsec] cross-track error
   double *terr,	// [arcsec] tangential-track error
   double *ragct,	// [deg] GCT RA from observer to object
   double *decgct);	// [deg] GCT Dec from observer to object


/* Return position and velocity with respect to the barycenter */
void puma_bary(
   void *p,		// puma handle returned by puma_init()
   int which,		// detection (1:ndet) to evaluate
   int integrate,       // [0/1] normal puma eval or true integration?  **NOTE!
   double *mjd,	        // MJD of this obs
   double dX[3],	// [m] object-barycenter displacement
   double dV[3],	// [m/s] object-barycenter velocity
   double px[3],        // [m,rad] r,theta,phi object-barycenter displacement
   double pv[3],        // [m/s,dad] vr,vtheta,vphi object-barycenter velocity
   double *sigma);      // [rad] point uncertainty
// **NOTE! 'DIFSYS=difsys.o' Makefile must enabled in order to use integrate=1

/* Return information from detection array for a particular obs/eph */
void puma_obsinfo(
   void *p,		// puma handle returned by puma_init()
   int which,		// detection (1:ndet) to evaluate
   int *eph,		// 0/1 for fitted data or ephemerides only
   double *t0,		// [day] origin MJD
   double *mjd,		// [day] MJD provided
   double *ra,		// [deg] RA provided (or RA_GCT for ephem)
   double *dec);	// [deg] Dec provided (or Dec_GCT for ephem)


/* Return selected Keplerian orbit parameters wrt Sun or Earth */
void puma_kepler(
   void *p,		// puma handle returned by puma_init()
   int earthctr,	// [0,1] wrt Sun (0) or Earth (1)
   double *E,		// [(km/s)^2] orbit energy
   double *ecc,		// [] orbit eccentricity
   double *incl,	// [deg] orbit inclination (ecliptic!)
   double *a,		// [AU,km] orbit semi-major axis
   double *period);	// [year,day] orbit period 

double puma_perigee(	// return distance of closest approach
   void *p,		// puma handle returned by puma_init()
   double *mjd,		// [day] MJD of closest approach
   double *ra,		// [deg] RA of closest approach
   double *dec,		// [deg] Dec of closest approach
   double *vtan);	// [deg/day] angular velocity at closest approach

/* Save or retrieve the obsid field from the puma->det->obs structure */
void puma_id(void *p, int idx, int put, char *id);
