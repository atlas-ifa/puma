/* Header file for orbit.c */

#define JT		/* Define to mess around! */

#define ABS(a) (((a) > 0) ? (a) : -(a))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define NINT(a)  (((a) > 0) ? (int)((a)+0.5) : (int)((a)-0.5))

/* Definition of REAL and associated math functions */
#define LONGREAL

/* 128 bit floating point */
#ifdef LONGREAL
typedef long double REAL;	/* 16 byte floating point */
#define SQRT(x) (sqrtl(x))
#define POW(x,y) (powl(x,y))
#define LOG(x) (logl(x))
#define EXP(x) (expl(x))
#define COS(x) (cosl(x))
#define SIN(x) (sinl(x))
#define TAN(x) (tanl(x))
#define ACOS(x) (acosl(x))
#define ASIN(x) (asinl(x))
#define ATAN(x) (atanl(x))
#define ATAN2(x,y) (atan2l(x,y))
#define SINH(x) (sinhl(x))
#define COSH(x) (coshl(x))
#define ATANH(x) (atanhl(x))

#define MINEPS 1e-12		/* Minimum error for difsys */


/* 64 bit floating point */
#else
typedef double REAL;	/* 16 byte floating point */
#define SQRT(x) (sqrt(x))
#define POW(x,y) (pow(x,y))
#define LOG(x) (log(x))
#define EXP(x) (exp(x))
#define COS(x) (cos(x))
#define SIN(x) (sin(x))
#define TAN(x) (tan(x))
#define ACOS(x) (acos(x))
#define ASIN(x) (asin(x))
#define ATAN(x) (atan(x))
#define ATAN2(x,y) (atan2(x,y))
#define SINH(x) (sinh(x))
#define COSH(x) (cosh(x))
#define ATANH(x) (atanh(x))

#define MINEPS 1e-9		/* Minimum error */

#endif



/* Fundamental constants */
#define cLIGHT (299792458.0)	/* [m/s] */
#define LN10 (2.302585092994)	/* ln(10) */
#define PI (3.14159265358979323846264338327950288419716939937510582) /* pi */
#define SECDAY (86400)		/* Seconds in a day */
#define GMSUN (1.327127e20)	/* [m^3/s^2] */
#define GMEARTH (3.986004418e14) /* [m^3/s^2] */
#define GMMOON (4.90390e12)	/* [m^3/s^2] */
#define GMJUPITER (1.26683e17)	/* [m^3/s^2] */
#define AU (1.49598E11) 	/* [m] */
#define DMOON (3.8844e8) 	/* [m] Lunar orbit semi-major axis */
#define RSUN (6.958e8) 		/* [m] Radius of the Sun */
#define YEAR (365.256366)	/* year [solar day] */
#define AGEO (42241095.7)	/* [m] aGEO=(GMEARTH*SECDAY^2/4/PI^2)^1/3 */

/* The WGS84 model of the Earth */
#define WGS84_a (6378137.0)		/* [m] */
#define WGS84_finv (298.257223563)	/* [] inverse flattening */
#define WGS84_GM (3.986004418e14)	/* [m^3/s^2] */
#define WGS84_J2 (1.081874e-3)		/* [] scaled oblateness */
#define WGS84_J3 (-2.532e-6)		/* [] scaled N-S asymmetry */
#define WGS84_omega (7.292115e-5)	/* [rad/sec] */

/* Phase space is 6 degrees of freedom */
#define NDIFEQ  6
/* Plenty of room for detections */
#define MAXPT 100000	/* Max number of detections */

/* A 3-vector */
typedef struct {
   REAL x;	/* [m] x position wrt inertial Earth */
   REAL y;	/* [m] y position wrt inertial Earth */
   REAL z;	/* [m] z position wrt inertial Earth */
} VEC;

/* Vector functions */
/* Dot product of two vectors: A.B */
#define DOT(A,B) ( (A).x*(B).x+(A).y*(B).y +(A).z*(B).z )

/* Cross product of two vectors: AxB */
#define XCRS(A,B) ( (A).y*(B).z-(A).z*(B).y )
#define YCRS(A,B) ( (A).z*(B).x-(A).x*(B).z )
#define ZCRS(A,B) ( (A).x*(B).y-(A).y*(B).x )

/* Triple product of three vectors: AxB.C = C.AxB = A.BxC = B.CxA */
#define TPROD(A,B,C) ( ((A).y*(B).z-(A).z*(B).y)*(C).x+((A).z*(B).x-(A).x*(B).z)*(C).y+((A).x*(B).y-(A).y*(B).x)*(C).z )

/* Observed properties of a detection */
typedef struct {
   int coosys;		/* Coordinate system of detection */
   char id[65];		/* observation ID */
   REAL mjd;		/* [UTC-MJD] centerpoint of observation */
   REAL ltt;		/* [sec] light travel time */
   REAL lng;		/* [rad] Observatory E longitude */
   REAL lat;		/* [rad] Observatory latitude */
   REAL alt;		/* [m] Observatory altitude */
   REAL ra;		/* [rad] observation RA */
   REAL dec;		/* [rad] observation Dec */
   REAL lambda;		/* [rad] observation lambda */
   REAL beta;		/* [rad] observation beta */
   REAL x;		/* [pixel] observation x (TSK) */
   REAL y;		/* [pixel] observation y (TSK) */
   REAL lst;		/* [rad] local sidereal time */
   REAL pnr[9];		/* PNR rotation matrix from terrestrial to inertial at MJD */
   REAL ut1utc;		/* [sec] (UT1-UTC) */
   REAL etime;		/* [sec] observation exposure time */
   REAL len;		/* [pix] major axis (streak length) */
   REAL psf;		/* [pix] minor axis PSF */
   REAL pa;		/* [deg] observation streak PA */
   REAL mag;		/* [AB] how bright? */
   REAL dmag;		/* [AB] magnitude uncertainty */
   REAL scale;		/* ["/pix] plate scale */
   int tpherr;		/* tphot error code */
} OBS;

/* Working detection */
typedef struct {
   OBS obs;		/* Raw detection quantities */
   REAL wgt;		/* [rad^-2] Uncertainty variance in RA/Dec */
   VEC O;		/* [m] Location of observatory wrt Sun */
   VEC U;		/* [m/s] Velocity of observatory wrt Sun */
   VEC X;		/* [ecl] Unit vector towards object from observatory */
   VEC E;		/* [m] Location of Earth's barycenter wrt Sun */
   VEC G;		/* [m] Location of Earth's geocenter wrt Sun */
   VEC M;		/* [m] Location of Moon's center wrt Sun */
   VEC W;		/* [m/s] Velocity of barycenter wrt Sun */
   VEC A;		/* [m/s2] Acceleration of barycenter wrt Sun */
   VEC XB;		/* [ecl] Unit vector towards object from barycenter */
   VEC OB;		/* [m] Location of observatory wrt barycenter */
   VEC P;		/* [m] Orbit prediction of object location wrt Sun */
   VEC V;		/* [m/s] Orbit prediction of object velocity wrt Sun */
   VEC XSAVE;		/* Unit vector towards object from observatory */
   int xsave_valid;	/* 0/1 if XSAVE has valid data */
   int chunk;		/* Sequence chunk */
   int use;		/* Use for fit? */
   REAL r;		/* [m] Object-observer distance */
   REAL rho;		/* [m] Object-barycentric distance */
   REAL rhoacc;		/* [m/s^2] Object-barycentric distance accel */
   REAL rog;		/* [m] Observer-barycenter distance */
   REAL omega;		/* [rad/sec] bary-obs phidot */
   REAL rorb;		/* [rad] Orbit RA */
   REAL dorb;		/* [rad] Orbit Dec */
   REAL rdot;		/* [rad/sec] Orbit RA angular velocity */
   REAL hdot;		/* [rad/sec] Orbit HA angular velocity */
   REAL ddot;		/* [rad/sec] Orbit Dec angular velocity */
   REAL tresid;		/* [rad] Angle between GCT and data along track */
   REAL xresid;		/* [rad] Angle between GCT data across track */
   REAL dpara;		/* [rad] Angle between orbit and data along track */
   REAL dperp;		/* [rad] Angle between orbit and data across track */
   REAL epara;		/* [rad] Error in angle between orbit and data along track */
   REAL eperp;		/* [rad] Error in angle between orbit and data across track */
   REAL dang;		/* [rad] Angle between orbit and observed */
   REAL dsig;		/* Sigma difference between orbit and observed */
} DET;

/* Orbit */
typedef struct {
   VEC X;		/* [m] location of object wrt sun (ecliptic) */
   VEC V;		/* [m/s] velocity wrt sun */
   REAL t0;		/* [UTC-MJD] reference time */
   VEC L;		/* [m^2/s] angular momentum */
   VEC Peri;		/* [m] periapse */
   REAL pra;		/* Pole RA */
   REAL pdec;		/* Pole Dec */
   REAL E;		/* [m^2/s^2] energy = v^2/2-GM/r*/
   REAL ecc;		/* Eccentricity = sqrt(1+2EL^2/GM^2) */
   REAL a;		/* [m] semi-major axis = GM/(2E) */
   REAL incl;		/* [rad] inclination = acos(Lz/L) = 90 - pole Dec */
   REAL Omega;		/* [rad] ascending node = atan(Lx,-Ly) = pole RA + 90 */
   REAL omega;		/* [rad] argument of periapse */
   REAL nu;		/* [rad] true anomaly: periapse to curposn */
   REAL period;		/* [sec] period = 2pi sqrt(a^3/GM) */
   REAL deriv[NDIFEQ*NDIFEQ];	/* Derivatives of kepler wrt state vector: */
   /* deriv[i+NDIFEQ*j] = dKepler_j / dstate_i */
   /* where Kepler_j = (incl, Omega, ecc, period, omega, nu) */
   /*    and state_i = (x, y, z, vx, vy, vz) */
} ORBIT;

/* PUMA data and fit structure */
typedef struct {
   int ndet;	/* Number of detections (including ephemeris requests) */
   int nfit;	/* Number of detections that are fitted */
   int id0;	/* Index of origin detection */
   int nalloc;	/* Allocation count of detection array */
   DET *det;	/* Detections */
   int ndate;	/* Order of acceleration polynomial */
   REAL R;	/* [m] current barycentric range at origin of fit */
   REAL VR;	/* [m/s] current barycentric dr/dt at origin of fit */
   VEC A[5];	/* polynomial coefficients of position fit */
   VEC g[2];	/* great circle pole and origin */
   REAL phi[4];	/* [rad/sec^2] angle and velocity around great circle */
} PUMA;

#define PUMA_NALLOC 256	/* Number of detections allocated at once for puma */


/* Noise model: floor +Q+ 1pix*dmag +Q+ trail*dmag) */
#define DPARA (0.05)	/* Floor uncertainty along track [pix] */
#define DPERP (0.025)	/* Floor uncertainty across track [pix] */
#define SECPIX (19.38)	/* Canon 5D3 and 135f/2 plate scale */

/* How many iteration cycles permitted mini?  (grossly excessive, we trust) */
#define MAXMINITER 100

/* Acceleration terms */
#define SOLARCONST (1361.0)	/* [W/m^2] Power from sun */

// #define DEBUG1

/* Prototypes */

/* Read an input file, return number of detections and data */
void read_input(char *infile, DET **det, int *npt);

/* Compute cartesian positions for observer and observation */
void cartesian(int npt, DET *d);

/* Local sidereal time (rad) at modified Julian date mjd at E long (rad) */
REAL lst(REAL mjd, REAL elng);

/* Quicksort program for REALs, sort index */
int qsort2(int n, REAL *x, int *idx);

/* Acceleration felt by object [m/s^2] */
int accel(REAL t, REAL *y, REAL *ydot, REAL *aux);

/* Update X,V from initial conditions at tnow to tgoal */
int timestep(REAL tnow, VEC *X, VEC *V, REAL tgoal/*MJD*/);

/* Integrate a system of differential equations */
int difsys(int func(REAL, REAL *, REAL *, REAL *), int n, REAL *h, 
	   REAL *x, REAL *y, REAL eps, REAL *s, REAL *aux, int *newh, int *err);

/* Return precession, nutation, and rotation matrix from mjd to 2000.0 */
void precnutrot(REAL mjd, REAL dUT1, REAL *PNR);

/* Moon position from center of Earth */
void xyzmoon(REAL mjd, VEC *moon);

/* Moon position from center of Earth */
void lp_xyzmoon(REAL mjd, VEC *moon);

/* Low precision moon position from center of Earth: ecliptic coordinates */
void lp_xyzmoon_ecliptic(REAL mjd, VEC *moon);

/* Low precision formulae for the sun, from Almanac p. C24 (1990) */
void xyzsun(REAL mjd, VEC *sun);

/* Low precision formulae for the sun, from Almanac p. C24 (1990) */
void xyzsun_lpmoon(REAL mjd, VEC *sun);

/* Calculate moon and sun params and locs for the range mjd1-mjd2, origin mjd0 */
void perturb(REAL mjd0, REAL mjd1, REAL mjd2);

/* Calculate delta times from three observations and an orbit pole */
int tripdt(REAL pra, REAL pdec, REAL GM, DET *d1, DET *d2, DET *d3,
	   VEC *X,	/* [m] 3D position at d1 */
	   VEC *V,	/* [m/s] 3D velocity at d1 */
	   VEC *P,	/* unit vector of pole */
	   VEC *PERI,	/* unit vector of periapse */
	   REAL *dt1,	/* [s] time between d1 and d2 */
	   REAL *dt2);	/* [s] time between d1 and d3 */

/* Compute Keplerian orbital constants from state vector */
void keparams(ORBIT *orb, REAL GM);

/* Compute Keplerian elements and derivatives from state vector */
void kepler(ORBIT *orb, REAL GM);

/* Compute the EGM2008 potential and gradient at r,lam,phi to degree lmax */
//void acc_egm08(double r, double lng, double lat, int lmax, int mmax,
//	       double *V, double *Vr, double *Vphi, double *Vlam);
void acc_egm08(REAL r, REAL lng, REAL lat, int lmax, int mmax,
	       REAL *V, REAL *Vr, REAL *Vphi, REAL *Vlam);

