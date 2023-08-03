/* Header file for pumalink */
/* v1.1 - 220228 John Tonry */

#define SECDAY (86400)          /* Seconds in a day */
#define AU (1.495978707E11)     /* [m] */

#define RMIN (0.001)    /* [AU] minimum R for a fit */


#define ABS(x) ((x)<0?(-(x)):(x))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define NINT(a)  (((a) > 0) ? (int)((a)+0.5) : (int)((a)-0.5))

#define NDIM 3    /* 3+3 = 6D phase space */
#define MAXPAR 6  /* possible linear parameters: 1,s,w,w^2,sw,w^2 */
#define NPAR 3    /* current linear parameters: 1,s,w... */

/* Global variables */
extern int PUMA_VERB;

/* One detection with trd information */
typedef struct {
   char *id;            /* observation ID */
   double mjd;          /* [UTC-MJD] midpoint of observation */
   double ra;           /* [deg] observation RA */
   double dec;          /* [deg] observation Dec */
   float xerr;          /* [arcsec] cross-track error */
   float terr;          /* [arcsec] trangential-track error */
   float lng;           /* [deg] Observatory E longitude */
   float lat;           /* [deg] Observatory latitude */
   float elev;          /* [m] Observatory elevation about geoid */
   int dgrp;            /* group number (or -1) */
   int idx;             /* index among all detections, sorted by ID */
} DET;

/* One set of detections with pumalink ephemeris information */
typedef struct {
   int ndet;                 /* Number of detections */
   DET **d;                  /* Detections */
   double mjd;               /* average MJD of constituent detections */
   double dt;                /* [sec] typical span between constituents */
   float S[2*NDIM];          /* Central estimate for S: x0,w0,x1,w1,x2,w2 */
   float A[NPAR*2*NDIM];     /* 1,s,w,... coeff: Asx0,Awx0, Asv0,Awv0, ...) */
   float rms[NDIM*2];        /* rms of for each S component */
   float cov[NDIM];          /* covariance between x,v for each S coordinate */
   float sigma;              /* [rad] detection uncertainty */
} VEC;

/* A detection list */
typedef struct {
   int dgrp;            /* detection group ID: disjoint set of detections */
   int sgrp;            /* subgroup ID: consistent state vector at ref time */
   int ndet;            /* Number of participating detections */
   int nalloc;          /* Allocation size of DET* array */
   DET **d;             /* Sorted(!) array of detection pointers */
} DETLIST;

/* A puma solution from an arbitrary list of detections */
typedef struct {
   DETLIST det;         /* Participating detections */
   double chin;         /* log(chi^2/N) from puma fit */
   double chin1;        /* log(chi^2/N) for constant position */
   double chin2;        /* log(chi^2/N) for two constant positions */
   double rfit;         /* [m] fitted barycentric distance */
   double wfit;         /* [s^-1] fitted rdot/r */
   double dlnr;         /* uncertainty in lnr */
   double dwr;          /* uncertainty in puma w=rdot/r */
   double rwcov;        /* r,w covariance from puma */
   double X[3];         /* unit vector at ephemeris time */
   double W[3];         /* [deg/day] tangential velocity */
} GRP;

/* A group comprised of one pair of vectors with pumalink fit and puma fit */
typedef struct {
   VEC *v1;             /* Vector 1 */
   VEC *v2;             /* Vector 2 */
   double chi;          /* log(min chi^2) from pumalink */
   double s0;           /* [AU^-1] s location of chi_min */
   double w0;           /* [km/s/AU] w location of chi_min */
   double ds;           /* [AU^-1] s,w covariance */
   double dw;           /* [km/s/AU] s,w covariance from pumalink */
   double swcov;        /* s,w covariance */
   GRP grp;             /* puma fit group */
} QUAD;

/* All the processing parameters */
typedef struct {
   double mjd;          /* common ephemeris MJD */
   double omega;        /* [deg/day] maximum ang vel for a vector */
   double dtmax;        /* [day] max interval between obs that might link */
   double mfrac;        /* [day] rounding fraction for reference MJD */
   int integrate;       /* [0/1] normal puma eval or true integration? */
   int split;           /* [0/1] split detections on reference time? */
   int nr;              /* steps in lnr for test grid */
   int nv;              /* steps in v for test grid */
   double rmin;         /* [AU] minimum radius of grid */
   double rmax;         /* [AU] maximum radius of grid */
   double vmin;         /* [km/s] minimum velocity of grid */
   double vmax;         /* [km/s] maximum velocity of grid */
   int ICnr;            /* steps in 1/r for initial condition grid */
   int ICnv;            /* steps in v for initial condition grid */
   double ICrmin;       /* [AU] minimum radius of IC grid */
   double ICrmax;       /* [AU] maximum radius of IC grid */
   double ICvmin;       /* [km/s] minimum velocity of IC grid */
   double ICvmax;       /* [km/s] maximum velocity of IC grid */
   double vtmax;        /* [km/s] maximum physical velocity from grid */
   double wdtmin;       /* minimum w*Dt permited from grid */
   double chieph;       /* Maximum chi^2/N from puma for linear grid */
   double dxmax;        /* [] max projected unit vector diff to link */
   double dwmax;        /* [rad/sec] max projected ang vel diff to link */
   double chimax;       /* Maximum chi^2 from pumalink that we will keep */
   double chinmax;      /* Maximum chi^2/N from puma that we will keep */
   double chiprod;      /* Maximum (chi^2 * chi^2/N) that we will keep */
   double subgrpxtol;   /* [rad] subgroup FoF unit vector tolerance */
   double subgrpwtol;   /* [deg/day] subgroup FoF ang vel tolerance */
   int ndet;            /* Total number of detections */
} LINKPAR;


/* Quicksorts */
void tsort_d(int n, double *x, int *idx);
void tsort_i(int n, int *x, int *idx);
void tsort_I(int n, long long int *x, int *idx);
void tsort_str(int n, char **x, int *idx);
void tsort_dp(int n, double *x, void* *idx);
void tsort_ip(int n, int *x, void* *idx);
void tsort_strp(int n, char **x, void* *idx);

/* Read an input file of trd data */
int read_trd(char *infile, DET **det, LINKPAR *linkpar);

/* Provide an intermediate, average MJD */
double avg_mjd(int ndet1, DET *det1, int ndet2, DET *det2);

/* Assign a unique index to all the detections */
int det_index(LINKPAR *linkpar, int ndet1, DET *det1, int ndet2, DET *det2);

/* Identify pairs of detections within ang vel w[deg/day] of each other */
int link_det(LINKPAR *linkpar, int n, DET *det, int epoch, VEC **vec);

/* Read an input list of detections and create vectors */
int read_vec(char *infile, LINKPAR *linkpar, int ndet, DET *det, int nvec, VEC **vec);

/* Solve nxn matrix eq Y = M X, return X for ny versions (mangle M,Y!) */
int linsolven(int n, int ny, double *Y, double *M, double *X);

/* Angle betweeen two points on the sphere, all [deg] */
double sphdist(double r1, double d1, double r2, double d2);

/* Parse a CSV string with real variables, return array and number */
int parse_eph(char *arg, int maxvar, double *var, char **id,
              double *lng, double *lat, double *elev);

/* Project one vector to the ephemeris time */
int pumephem(LINKPAR *linkpar, void *puma, VEC *vec, FILE *fpscat, int nscat);

/* Return chi^2 from a pumalink assessment of a pair of vectors */
double pumalink(VEC *v1, VEC *v2, double ephmjd, double *sw);

/* Return chin from a puma fit to all four detections of a pair of vectors */
double pumatest(void *puma, int ndet, DET **d, LINKPAR *linkpar,
                double *X, double *Wtan, double *rv, double *chin1);

/* Match all vectors according to tol, test pumalink and puma quad */
int pumaquad(void *puma, int nvec, VEC *vec, LINKPAR *linkpar, QUAD **quad);

/* Group all detections that are linked by a quad */
int detgroup(int nquad, QUAD *quad);

/* Split each detection group into FoF subgroups with overlapping quads */
int quadgroup(int nquad, QUAD *quad, LINKPAR *linkpar);

/* Do a puma fit to all the detections of quad subgroups */
int pumamult(void *puma, int nquad, QUAD *quad, LINKPAR *linkpar, GRP **mult);
