/* Link trd data from two epochs */
/* v1.0 - 221226 John Tonry */
/* Syntax: pumalink trd [options] */
/* Subversion r.15472 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "puma.h"
#include "pumalink.h"
#include "kdtree.h"

static int VERBOSE=1;

#define NEPHDETMIN 100  /* vector det count to fit for pumephem grid */

#define DETPROX 0.003   /* [deg] unacceptable detection separation */

#define MAXCOV 0.999    /* Clip fit RMS covariance at MAXCOV */
#define MAXCRV 0.99999  /* Clip propagated r,v covariance at MAXCRV */
#define VARBLOAT 1.0    /* Bloat factor for obs sigma^2 for pumalink() */

extern int NXROCK;      /* libpuma count of calls to xrock() */

/* Heavy duty debugging */
char *DEBUG_FNAME="/tmp/pumalink_debug";
FILE *TRACK_SINGLE_DET=NULL;    /* Dump all detections to DEBUG_FNAME.alldet */
FILE *TRACK_SINGLE_PAIR=NULL;   /* Dump all pairs to DEBUG_FNAME.allpair */
FILE *TRACK_SINGLE_VEC=NULL;    /* Dump all vectors to DEBUG_FNAME.allvec */
FILE *TRACK_SINGLE_QUAD=NULL;   /* Dump all quads to DEBUG_FNAME.allquad */
FILE *DEBUG_QUAD=NULL;          /* Save all quads */

/* log10(x) but -9.999 for x<=0 */
#define LOGCLIP(x) ( ((x)>0) ? log((x))/log(10) : -9.999 )

void usage(int error) 
{
/* Emit usage message */
   fprintf(stderr, "Usage: pumavec TRD_file [options]\n");
   exit(error);
}

/* Execute to get a time of day */
#define TOD(msg) gettimeofday(&tv2, NULL); dt0 = (tv2.tv_sec - tv0.tv_sec) + 1e-6*(tv2.tv_usec - tv0.tv_usec); dt1 = (tv2.tv_sec - tv1.tv_sec) + 1e-6*(tv2.tv_usec - tv1.tv_usec); fprintf(stderr, "t= %8.3f  dt= %8.3f %s\n", dt0, dt1, msg); tv1 = tv2;


int main(int argc, char **argv)
{
   int i, j, ndet, nvec, nscat, nquad, ngrp, nmult, maxsubgrp;
   int nfailpt, nfailfit, DEBUG;
   double dt0, dt1, trmax;

   char *trd;           /* Input detection files */
   char *vf1, *vf2;     /* Input vector specification files */
   char *quadout, *multout;     /* Output files for quads and groups */
   char *scatout, *dbquad;      /* Output files for diagnostics */
   FILE *fpscat, *fpquad, *fpmult;
   void *puma;
   DET *det;            /* Detection array */
   VEC *vec;            /* Vector array */
   int nvec1, nvec2;
   QUAD *quad;          /* Quad array */
   GRP *mult;          /* Group array */
   LINKPAR linkpar;     /* link parameters */

   char info[1024];
   struct timeval tv0, tv1, tv2;

/* Start time */
   gettimeofday(&tv0, NULL);
   tv1 = tv0;

/* Parse the arguments */
   if (argc < 2) usage(1);
   trd = argv[1];

/* Defaults */
   vf1 = NULL;          /* vector input lists for epoch 1 */
   vf2 = NULL;          /* vector input lists for epoch 2 */
//   vec1 = vec2 = NULL;

   linkpar.omega = 5;      /* [deg/day] maximum ang vel for a vector */
   linkpar.dtmax = 0.1;    /* [day] max interval between obs that might link */
   linkpar.mfrac = 0.001;  /* [day] rounding fraction for refernce MJD */
   linkpar.split = 1;      /* split detections on reference time? */
   linkpar.mjd = 0;        /* [MJD] split date */
   linkpar.integrate = 0;  /* [0/1] normal puma eval or true integration */

/* Nominal, modest grid for fit */
   linkpar.nr = 5;           /* Number of grid search in r */
   linkpar.nv = 5;           /* Number of grid search in vr */
   linkpar.rmin = 0.02;      /* [AU] Minimum r */
   linkpar.rmax = 4;         /* [AU] Maximum r */
   linkpar.vmin = -20;       /* [km/s] Minimum radial vr */
   linkpar.vmax = 20;        /* [km/s] Maximum radial vr */
   trmax = 10;             /* [AU] triangle grid rmax */

   linkpar.wdtmin = 0.5;     /* [] minimum w*Dt for linear fit input */
   linkpar.vtmax = 100;      /* [km/s] max V permitted for linear fit input */
   linkpar.chieph = 100;     /* Maximum chi^2/N for grid test points */

// Use a triangular s,w grid by default?
#if 0
   linkpar.nr = 5;           /* Number of grid search in r */
   linkpar.nv = 0;           /* Number of grid search in vr */
   linkpar.rmin = 0.1;       /* [AU] Minimum r */
   linkpar.rmax = 4;         /* [AU] Maximum r */
   linkpar.vmin = -20;       /* [km/s] Minimum radial vr */
   linkpar.vmax = 20;        /* [km/s] Maximum radial vr */
#endif

/* Nominal, modest grid for fit */
   linkpar.ICnr = 1;         /* Number of IC grid search in r */
   linkpar.ICnv = 1;         /* Number of IC grid search in vr */
   linkpar.ICrmin = 0.36;    /* [AU] Minimum r for IC */
   linkpar.ICrmax = 0.36;    /* [AU] Maximum r for IC */
   linkpar.ICvmin = 0;       /* [km/s] Minimum radial vr for IC */
   linkpar.ICvmax = 0;       /* [km/s] Minimum radial vr for IC*/

   linkpar.dxmax = 0.004;    /* [] max projected unit vector diff to link */
   linkpar.dwmax = 1e-7;     /* [rad/sec] max projected ang vel diff to link */

//   linkpar.chimax =  200; /* Maximum chi^2 from pumalink that we will keep */
   linkpar.chimax =  25;    /* Maximum chi^2 from pumalink that we will keep */
//   linkpar.chinmax =  30;  /* Maximum chi^2/N from puma that we will keep */
   linkpar.chinmax =  25;    /* Maximum chi^2/N from puma that we will keep */
   linkpar.chiprod = 200;    /* Maximum (chi^2 * chi^2/N) that we will keep */

   linkpar.subgrpxtol = 0.0002;  /* [rad] subgroup FoF unit vector tolerance */
   linkpar.subgrpwtol = 0.002;   /* [deg/day] subgroup FoF ang vel tolerance */

/* Output defaults */
   quadout = NULL;	/* quad output file */
   fpquad = stdout;
   multout = NULL;	/* multiple output file */
   fpmult = stdout;

/* Diagnostic puma output to a scatter file? */
   fpscat = NULL;	/* Statistics scatter output file */
   scatout = NULL;	/* Statistics scatter output file */
   nscat = 0;           /* mess with the data this many times */
   dbquad = NULL;	/* quad debug file */

/* Deep debugging? */
   DEBUG = 0;

/* Parse the arguments */
   for(i=2; i<argc; i++) {
      if (strcmp(argv[i], "-verbose") == 0) {	   /* debug? */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%d", &PUMA_VERB);
	 VERBOSE = PUMA_VERB;

      } else if (strcmp(argv[i], "-DEBUG") == 0) {	   /* debug? */
         DEBUG = 1;

      } else if (strcmp(argv[i], "-vec1") == 0) {  /* Vector lists epoch 1 */
         if (i == argc - 1) usage(1);
	 vf1 = argv[++i];

      } else if (strcmp(argv[i], "-vec2") == 0) {  /* Vector lists epoch 1 */
         if (i == argc - 1) usage(1);
	 vf2 = argv[++i];

      } else if (strcmp(argv[i], "-omega") == 0) {/* [deg/day] max ang vel */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &linkpar.omega);

      } else if (strcmp(argv[i], "-dtmax") == 0) {/* [day] max dt for a pair */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &linkpar.dtmax);

      } else if (strcmp(argv[i], "-mfrac") == 0) {/* [day] rounding for teph */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &linkpar.mfrac);

      } else if (strcmp(argv[i], "-split") == 0) {/* split detections? */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%d", &linkpar.split);

      } else if (strcmp(argv[i], "-mjd") == 0) {  /* split MJD? */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &linkpar.mjd);

      } else if (strcmp(argv[i], "-integrate") == 0) {  /* puma integrate? */
         linkpar.integrate = 1;

      } else if (strcmp(argv[i], "-dxmax") == 0) {/* [] unit vector diff */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &linkpar.dxmax);

      } else if (strcmp(argv[i], "-dwmax") == 0) {/* [rad/s] ang vel diff */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &linkpar.dwmax);

      } else if (strcmp(argv[i], "-chimax") == 0) {/* maximum chi^2 tolerated */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &linkpar.chimax);

      } else if (strcmp(argv[i], "-chinmax") == 0) {/* maximum chi^2/N */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &linkpar.chinmax);

      } else if (strcmp(argv[i], "-chiprod") == 0) {/* maximum chi*chin */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &linkpar.chiprod);

      } else if (strcmp(argv[i], "-wdt") == 0) {/* maximum wDt tolerated */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &linkpar.wdtmin);

      } else if (strcmp(argv[i], "-vmax") == 0) {/* [km/s] max V tolerated */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &linkpar.vtmax);

      } else if (strcmp(argv[i], "-grid") == 0) {  /* grid spec */
	 linkpar.nr = linkpar.nv = 20;  /* Grid count for nominal limits above */
	 if(i<argc-1 && sscanf(argv[i+1], "%d,%lf,%lf,%d,%lf,%lf",
                &linkpar.nr, &linkpar.rmin, &linkpar.rmax, &linkpar.nv, &linkpar.vmin, &linkpar.vmax) == 6) {
	    i++;
	 }

      } else if (strcmp(argv[i], "-chieph") == 0) {/* max chi^2/N for grid */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &linkpar.chieph);

      } else if (strcmp(argv[i], "-trmax") == 0) {/* [AU] rmax for trigrid */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf", &trmax);

      } else if (strcmp(argv[i], "-trig") == 0) {  /* triangualar grid spec */
	 linkpar.nr = linkpar.nv = 20;  /* Grid count for nominal limits above */
	 if(i<argc-1 && sscanf(argv[i+1], "%d,%lf,%lf",
                &linkpar.nr, &linkpar.rmin, &linkpar.vmax) == 3) {
            linkpar.nv = 0;
            linkpar.rmax = trmax;
            linkpar.vmin = -linkpar.vmax;
	    i++;
	 }

      } else if (strcmp(argv[i], "-grptol") == 0) {/* [rad,deg/day] subgroup tol */
         if (i == argc - 1) usage(1);
         sscanf(argv[++i], "%lf,%lf", &linkpar.subgrpxtol, &linkpar.subgrpwtol);

      } else if (strcmp(argv[i], "-quad") == 0) {  /* Quad out file */
         if (i == argc - 1) usage(1);
	 quadout = argv[++i];

      } else if (strcmp(argv[i], "-grp") == 0) {  /* Group out file */
         if (i == argc - 1) usage(1);
	 multout = argv[++i];

      } else if (strcmp(argv[i], "-nscat") == 0) {  /* [km/s] linear initial vr */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%d", &nscat);

      } else if (strcmp(argv[i], "-scat") == 0) {  /* Grid scatter out file */
         if (i == argc - 1) usage(1);
	 scatout = argv[++i];

      } else if (strcmp(argv[i], "-dbquad") == 0) {  /* DEBUG quads */
         if (i == argc - 1) usage(1);
	 dbquad = argv[++i];

      } else {
	 fprintf(stderr, "Unknown arg `%s'\n", argv[i]);
	 exit(1);
      }
   }

/* Open output quad file if requested */
   if(quadout != NULL) {
      if(strcmp(quadout, "-") == 0) fpquad = stdout;
      else if((fpquad = fopen(quadout, "w")) == NULL) {
	 fprintf(stderr, "Cannot open %s for writing quads\n", quadout);
	 exit(1);
      }
   }

/* Open output file for multiple detection groups if requested */
   if(multout != NULL) {
      if(strcmp(multout, "-") == 0) fpmult = stdout;
      else if((fpmult = fopen(multout, "w")) == NULL) {
	 fprintf(stderr, "Cannot open %s for writing groups\n", multout);
	 exit(1);
      }
   }

/* Open output ephemeris scatter file if requested */
   if(scatout != NULL) {
      if(strcmp(scatout, "-") == 0) fpscat = stdout;
      else if((fpscat = fopen(scatout, "w")) == NULL) {
	 fprintf(stderr, "Cannot open %s for writing\n", scatout);
	 exit(1);
      }
   }

/* Open quad debug file if requested */
   if(dbquad != NULL && (DEBUG_QUAD = fopen(dbquad, "w")) == NULL) {
      fprintf(stderr, "Cannot open %s for writing\n", dbquad);
      exit(1);
   }
   
/* Open deep debug files if requested */
   if(DEBUG) {
      sprintf(info, "%s.alldet", DEBUG_FNAME);
      TRACK_SINGLE_DET=fopen(info, "w");
      sprintf(info, "%s.allpair", DEBUG_FNAME);
      TRACK_SINGLE_PAIR=fopen(info, "w");
      sprintf(info, "%s.allvec", DEBUG_FNAME);
      TRACK_SINGLE_VEC=fopen(info, "w");
      sprintf(info, "%s.allquad", DEBUG_FNAME);
      TRACK_SINGLE_QUAD=fopen(info, "w");
   }



/* Read the trd file */
   ndet = read_trd(trd, &det, &linkpar);
   sprintf(info, "  Read %d detections from %s, <mjd> %.3f", ndet, trd, linkpar.mjd);
   TOD(info);

/* Link detections if both epochs do not have a vector list provided */
   nvec = 0;
   vec = NULL;
   if(vf1 == NULL || vf2 == NULL) {
      i = 0;
      if(vf2 != NULL) i = 1;
      if(vf1 != NULL) i = 2;
      nvec = link_det(&linkpar, ndet, det, i, &vec);
      sprintf(info, "  Detections form %d vectors within %.1f deg/day and %.2f day",
              nvec, linkpar.omega, linkpar.dtmax);
      TOD(info);
   }

/* Add vectors provided for epoch 1 */
   if(vf1 != NULL) {
      nvec1 = read_vec(vf1, &linkpar, ndet, det, nvec, &vec);
      nvec += nvec1;
      sprintf(info, "  User first epoch %d vectors from %s", nvec1, vf1);
      TOD(info);
   }

/* Add vectors provided for epoch 2 */
   if(vf2 != NULL) {
      nvec2 = read_vec(vf2, &linkpar, ndet, det, nvec, &vec);
      nvec += nvec2;
      sprintf(info, "  User second epoch %d vectors from %s", nvec2, vf2);
      TOD(info);
   }

/* Initialize puma */
   if(puma_init(&puma)) {
      fprintf(stderr, "Cannot initialize puma\n");
      exit(1);
   }

/* Project vectors to reference time */
   nfailpt = nfailfit = 0;
   for(i=0; i<nvec; i++) {
      j = pumephem(&linkpar, puma, vec+i, fpscat, nscat);
      if(j == -1) nfailfit++;
      if(j == -2) nfailpt++;
      if((i%10000) == 0) fprintf(stderr, "eph %8d %8d %8d  \r",
                                 i, nfailpt, nfailfit);
   }
   sprintf(info, "  %d vectors extrapolated with %d extrap, %d fit failures",
           nvec, nfailpt, nfailfit);
   TOD(info)

/* If a scatter file is requested, quit at this point */
   if(fpscat != NULL) {
      sprintf(info, "\nScatter file `%s' written, now exit...", scatout);
      TOD(info);
      return(0);
   }

/* Test all quads */
   nquad = pumaquad(puma, nvec, vec, &linkpar, &quad);
   sprintf(info, "  Vectors link to %d quads", nquad);
   TOD(info);

/* Group all quads that share a detection */
   ngrp = detgroup(nquad, quad);
   
   sprintf(info, "  Quads group into %d detection groups", ngrp);
   TOD(info);

/* Subgroup all detection groups by unit vector */
   maxsubgrp = quadgroup(nquad, quad, &linkpar);

   sprintf(info, "  Quad state vector subgroup max size %d", maxsubgrp);
   TOD(info);

/* Do a puma fit to all the detections of quad subgroups */
   nmult = pumamult(puma, nquad, quad, &linkpar, &mult);
   sprintf(info, "  Quads group into %d state vector subgroups", nmult);
   TOD(info);


/* Put out all quads, sorted by detection group */
   QUAD **squad = (QUAD **)calloc(nquad, sizeof(QUAD *));
   double *srtgrp = (double *)calloc(nquad, sizeof(double));

   for(i=0; i<nquad; i++) {
      squad[i] = (void *)(quad+i);
      srtgrp[i] = quad[i].grp.det.dgrp + quad[i].grp.det.sgrp/(double)(maxsubgrp+1);
   }
   tsort_dp(nquad, srtgrp, (void **)squad);

/* Write all the quads */
   int dgrp=0, sgrp=0, qgrp;
   char *grpfmt="%4d.%-4d %6.3f %6.3f %6.3f %6.3f %7.4f %8.5f %5.3f %5.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %3d ";
   char *snglfmt="%4d.0%-3d %6.3f %6.3f %6.3f %6.3f %7.4f %8.5f %5.3f %5.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %3d ";
   char *fmt;

   for(i=0; i<nquad; i++) {
/* Reset detection group and singleton group counter */
      if(squad[i]->grp.det.dgrp != dgrp) {
         dgrp = squad[i]->grp.det.dgrp;
         sgrp = 0;
      }
/* Singleton group gets a state group with a leading zero */
      if(squad[i]->grp.det.sgrp == 0) {
         sgrp += 1;
         qgrp = sgrp;
         fmt = snglfmt;
/* Mult group gets a state group with no leading zero */
      } else {
         qgrp = squad[i]->grp.det.sgrp;
         fmt = grpfmt;
      }

/* For singleton quads construct a running singleton ID */
      fprintf(fpquad, fmt, squad[i]->grp.det.dgrp, qgrp,
              squad[i]->chi, squad[i]->grp.chin,
              squad[i]->grp.chin1, squad[i]->grp.chin2,
              squad[i]->grp.rfit/AU, squad[i]->grp.wfit*SECDAY,
              squad[i]->grp.dlnr, squad[i]->grp.dwr,
              squad[i]->grp.X[0], squad[i]->grp.X[1], squad[i]->grp.X[2],
              squad[i]->grp.W[0], squad[i]->grp.W[1], squad[i]->grp.W[2],
              squad[i]->grp.det.ndet);
      for(j=0; j<squad[i]->grp.det.ndet; j++) fprintf(fpquad, "%s%c",
                        squad[i]->grp.det.d[j]->id,
                        (j==squad[i]->grp.det.ndet-1) ? '\n' : ',');
   }

/* Put out all quad groups */
   for(i=0; i<nmult; i++) {
      fprintf(fpmult, "%4d.%-4d %6.3f %6.3f %6.3f %6.3f %7.4f %8.5f %5.3f %5.3f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %3d ",
              mult[i].det.dgrp, mult[i].det.sgrp,
              0.0, mult[i].chin, mult[i].chin1, mult[i].chin2,
              mult[i].rfit/AU, mult[i].wfit*SECDAY, mult[i].dlnr, mult[i].dwr,
              mult[i].X[0], mult[i].X[1], mult[i].X[2],
              mult[i].W[0], mult[i].W[1], mult[i].W[2], mult[i].det.ndet);
      for(j=0; j<mult[i].det.ndet; j++) fprintf(fpmult, "%s%c",
                        mult[i].det.d[j]->id, (j==mult[i].det.ndet-1) ? '\n' : ',');
   }


   exit(0);
}

/* Do a small grid search then a puma fit */
double puma_gridfit(void *puma,  double *lnr, double *wr, double *dlnr, double *dwr, double *rwcov, LINKPAR *linkpar)
{
   int i, k;
   double chin=1e9, tchin, tlnr, twr, trms[4];

/* Grid search for best lnr, wr using linkpar */
   for(k=0; k<linkpar->ICnv; k++) {
      for(i=0; i<linkpar->ICnr; i++) {
         tlnr = log(linkpar->ICrmin) +                 // [lnAU]
            log(linkpar->ICrmax/linkpar->ICrmin) *
            i / (double)(MAX(1,linkpar->ICnr-1));
         twr = (linkpar->ICvmin +                      // [/day]
                (linkpar->ICvmax-linkpar->ICvmin) *
                k / (double)(MAX(1,linkpar->ICnv-1))) *
               1e3 / (exp(tlnr)*AU) * SECDAY;
         tchin = puma_chi(puma, tlnr, twr, trms);
         if((i==0 && k==0) || tchin < chin) {
            chin = tchin;
            *lnr = tlnr;
            *wr = twr;
         }
      }
   }

/* Fit for lnr, wr using the best lnr,wr */
   chin = puma_fit(puma, lnr, wr, dlnr, dwr, rwcov);
   return(chin);
}


/* Fit the detections from one array of quad pointers */
int multfit(void *puma, int n, QUAD **q, LINKPAR *linkpar, GRP *m)
{
   int i, k, ndet, pumid, nn;
   double lnr, wr, dlnr, dwr, rwcov, sigma, rho, rhodot, mjd, ephmjd;
   double X[3], V[3], dX[3], dV[3], dr=atan(1)/45;
   int *idx;
   DET **d;

   ephmjd = linkpar->mjd;
   
/* How many potential detections? */
   for(k=ndet=0; k<n; k++) ndet += q[k]->grp.det.ndet;

   d = (DET **)calloc(ndet, sizeof(DET *));
   idx = (int *)calloc(ndet, sizeof(int));

/* Load up the IDs and detections, sort, then strip the dups */
   for(k=nn=0; k<n; k++) {
      for(i=0; i<q[k]->grp.det.ndet; i++) {
         d[nn] = q[k]->grp.det.d[i];
         idx[nn++] = q[k]->grp.det.d[i]->idx;
      }
   }
/* Sort on detection index (MJD,ID) */
   tsort_ip(nn, idx, (void **)d);
   for(i=ndet=1; i<nn; i++) if(idx[i] != idx[i-1]) d[ndet++] = d[i];

/* Save the number of detections, create detection pointer array */
   m->det.ndet = ndet;
   m->det.d = (DET **)calloc(ndet, sizeof(DET *));
   m->det.dgrp = q[0]->grp.det.dgrp;
   m->det.sgrp = q[0]->grp.det.sgrp;

/* Do a puma fit */
   puma_reset(puma);

/* Load the detections into puma */
   for(i=0; i<ndet; i++) {
      m->det.d[i] = d[i];
      pumid = puma_obs(puma, 0, d[i]->mjd, d[i]->ra, d[i]->dec,
                   (double)d[i]->xerr, (double)d[i]->terr,
                   (double)d[i]->lng, (double)d[i]->lat, (double)d[i]->elev);
   }

/* Fit for lnr, wr with grid search */
   m->chin = puma_gridfit(puma, &lnr, &wr, &dlnr, &dwr, &rwcov, linkpar);

/* Test single or double constant position detections */
   double chin2[2];
   m->chin1 = puma_const(puma, chin2);
   m->chin2 = chin2[0];

/* Convert chi^2 to log10 and clip */
   m->chin = LOGCLIP(m->chin);
   m->chin1 = LOGCLIP(m->chin1);
   m->chin2 = LOGCLIP(m->chin2);

/* Evaluate barycenter location at ephmjd */
   pumid = puma_obs(puma, 1, ephmjd, 0., 0., 0., 0., 0., 0., 0.);
   puma_bary(puma, pumid, linkpar->integrate, &mjd, X, V, dX, dV, &sigma);

/* At the requested epoch... */
   rho = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
   rhodot = (X[0]*V[0]+X[1]*V[1]+X[2]*V[2]) / rho;

/* Unit vector */
   X[0] /= rho;
   X[1] /= rho;
   X[2] /= rho;

/* Return particulars of the fit */
   m->rfit = rho;               // [m]
   m->wfit = rhodot/rho;        // [sec^-1]
   m->dlnr = dlnr;
   m->dwr = dwr;
   m->rwcov = rwcov;
   
   m->X[0] = X[0];
   m->X[1] = X[1];
   m->X[2] = X[2];

/* Vtan = tangential angular velocity [rad/sec] */
   m->W[0] = ((V[0] - X[0] * rhodot) / rho) / dr * SECDAY;
   m->W[1] = ((V[1] - X[1] * rhodot) / rho) / dr * SECDAY;
   m->W[2] = ((V[2] - X[2] * rhodot) / rho) / dr * SECDAY;

   free(idx);
   free(d);

   return(ndet);
}


#define GRPALLOC 10

/* Do a puma fit to all the detections of quad subgroups */
int pumamult(void *puma, int nquad, QUAD *quad, LINKPAR *linkpar, GRP **mult)
{
   int nmult, nalloc, n, k, nqok, ndet=0;
   double maxsub;

/* List of pointers to all quads */
   QUAD **q = (QUAD **)calloc(nquad, sizeof(QUAD *));
/* Groups numbers of all quads */
   double *grp = (double *)calloc(nquad, sizeof(double));

   nmult = nalloc = 0;
   *mult = NULL;

/* Build a group index and quad pointers */
   for(k=0, maxsub=0.0; k<nquad; k++)
      if(quad[k].grp.det.sgrp > maxsub) maxsub = quad[k].grp.det.sgrp;
   maxsub = exp(log(10)*ceil(log(maxsub+1)/log(10)));  // next power of 10
   if(maxsub < 1) maxsub = 1;

/* Get the lists of groups, leaving out bad or ungrouped quads */
   for(k=nqok=0; k<nquad; k++) {
/* Inhibit singleton and bad quads */
      if(quad[k].grp.det.dgrp > 0 && quad[k].grp.det.sgrp > 0) {
         grp[nqok] = quad[k].grp.det.dgrp + quad[k].grp.det.sgrp / maxsub;
         q[nqok] = quad + k;
         nqok++;
      }
   }
/* Sort all the quads by group number */
   tsort_dp(nqok, grp, (void **)q);

/* Collect each quad subgroup, make a GRP, and fit its list of detections*/
   for(k=0; k<nqok; k+=n) {
      for(n=1; n<nqok-k && grp[k+n] == grp[k+n-1]; n++);
      if(n == 1) continue;

/* Ensure there is space for the result */
      if(nmult >= nalloc) {
         nalloc += GRPALLOC;
         *mult = (GRP *)realloc(*mult, nalloc*sizeof(GRP));
      }

      ndet += multfit(puma, n, q+k, linkpar, (*mult)+nmult);
                     
      nmult++;
   }

   free(q);
   free(grp);

   return(nmult);
}


#define DGALLOC 10 

/* Insert the sorted entries of d1 into sorted d2, return number of unique */
void detlist_insert(DETLIST *d1, DETLIST *d2)
{
   int i, i1, i2;

/* Ensure there is space to add all these detections to the list */
   if(d2->ndet+d1->ndet >=  d2->nalloc) {
      d2->nalloc = d2->ndet+d1->ndet + DGALLOC;
      d2->d = (DET **)realloc(d2->d, d2->nalloc*sizeof(DET *));
   }

/* Insert detection pointers from d1 into d2, maintaining sort */
   for(i1=i2=0; i1<d1->ndet; i1++) {
      for(  ; (i2 < d2->ndet) && (d2->d[i2]->idx < d1->d[i1]->idx); i2++);
      if(i2 < d2->ndet) {
         if(d1->d[i1]->idx == d2->d[i2]->idx) continue;
         for(i=d2->ndet; i>=i2; i--) d2->d[i+1] = d2->d[i];
      }
      d2->d[i2] = d1->d[i1];
      d2->ndet += 1;
   }
   return;
}

/* Test all the detections of q against all the detections of dg
 * If q and dg have detections from from the same observation (same MJD)
 * which are excessively far away from one another to plausibly be the
 * same object, return 0: i.e. not consistent_det.
 * If there is no inconsistent detection from the same observation,
 * sort all the new detections from q into the list of detections of dg,
 * and return 1: i.e. consistent_det.
 */

/* Are the detections of this quad consistent with the extant list? */
int consistent_det(QUAD *q, DETLIST *dg)
{
// NO PRETENSE of EFFICIENCY!
   int i, j;
   DET *d;
   double dth;

/* Summed distance between each detection from the quad and the list */
   for(i=0, dth=0.0; i<dg->ndet; i++) {
      for(j=0; j<q->grp.det.ndet; j++) {
         d = q->grp.det.d[j];
/* Same observation?  If not it cannot contribute to the distance */
         if(d->mjd != dg->d[i]->mjd) continue;  // not same obs: OK
         dth += sphdist(d->ra, d->dec, dg->d[i]->ra, dg->d[i]->dec);
      }
   }

// FIXME: should perhaps test for different detection methods on same obs
// FIXME: DETPROX should be exposed to the user

/* Too far: cannot be consistent */
   if(dth > DETPROX) return(0);

/* Insert these detections into the detection list */
   detlist_insert(&(q->grp.det), dg);

   return(1);
}

// #define SUBGRPTEST   /* lots of debug output during consolidate_subgrps() */


/* grp is an array of ngrp DETLISTs that list the detections for
 * each subgroup of a detection group.  q is an array of nquad quad
 * pointers that belong in the same detection group; the q with
 * subgroup 0 were never matched to a subgroup within reference
 * location tolerance.
 *
 * The function of consolidate_subgrp() is identify all subgroups
 * whose detections lie entirely within some other subgroup, and wipe
 * out that subgroup, assigning the larger group's subgroup number to
 * the quad's of the smaller subgroup.  Singleton quads (one with no
 * subgroup assignment) which did not spatially match into a group are
 * also considered valid subgroups for dissolution.
 */

/* Rename the quads of any subgroup that is fully contained within another */
void consolidate_subgrps(LINKPAR *linkpar, int ngrp, DETLIST *grp,
                         int nquad, QUAD **qptr)
{
   int i, j, k, n, l, m, idet;
   DETLIST *dptr;

/* Pointers to all the groups, sorted biggest to smallest */
   DETLIST **grpdet = (DETLIST **)calloc(ngrp+nquad, sizeof(DETLIST *));
   int *gsize = (int *)calloc(ngrp+nquad, sizeof(int));

/* Sort pointers to the groups and quads from largest to smallest group size */
   for(i=n=0; i<ngrp; i++, n++) {
      grpdet[n] = &grp[i];
      gsize[n] = -grp[i].ndet;
   }
   for(i=0; i<nquad; i++) {
/* Skip groups of size 4 since they cannot contain any other group */
      if(qptr[i]->grp.det.ndet <= 4) continue;
      grpdet[n] = &qptr[i]->grp.det;
      gsize[n++] = -qptr[i]->grp.det.ndet;
   }
/* Sort the grpdet array from biggest to smallest DETLIST */
   tsort_ip(n, gsize, (void **)grpdet);

/* Make an array for all the detections and grp/quad that claim them */
/* Depends on 0 initialization by calloc()! */
   int *ngpd = (int *)calloc(linkpar->ndet, sizeof(int)); // (n group per det)
   int *nalloc = (int *)calloc(linkpar->ndet, sizeof(int));
   DETLIST ***detgrp = (DETLIST ***)calloc(linkpar->ndet, sizeof(DETLIST **));

/* Load the detection array with all the groups or quads that claim them */
   for(k=0; k<ngrp+nquad; k++) {
      if(k < ngrp) dptr = grp + k;
      else         dptr = &(qptr[k-ngrp]->grp.det);
      for(i=0; i<dptr->ndet; i++) {
         idet = dptr->d[i]->idx;
         if(ngpd[idet] >= nalloc[idet]) {
            nalloc[idet] = ngpd[idet] + 10;
            detgrp[idet] = (DETLIST **)realloc(detgrp[idet],
                                               nalloc[idet]*sizeof(DETLIST *));
         }
         detgrp[idet][ngpd[idet]] = dptr;
         ngpd[idet] += 1;
      }
   }

/* Visit each group biggest to smallest: grpdet[k] */
   for(k=0; k<n; k++) {

/* For each detection in this group: grpdet[k]->d[j] */
      for(j=0; j<grpdet[k]->ndet; j++) {

/* Visit each group that shares this detection */
         idet = grpdet[k]->d[j]->idx;
         for(i=0; i<ngpd[idet]; i++) {

/* Test d=detgrp[idet][i] for being a subset of grpdet[k] */
            dptr = detgrp[idet][i];
/* Skip if too big or the same group */
            if(dptr->ndet > grpdet[k]->ndet || dptr == grpdet[k]) continue;

/* Test all detections of the smaller group against the bigger */
            for(l=0; l<dptr->ndet; l++) {  // not nec sorted: consistent_det()
               for(m=0; m<grpdet[k]->ndet && 
                        dptr->d[l] != grpdet[k]->d[m]; m++);
               if(m == grpdet[k]->ndet) break; // matched no det in grpdet[k]
            }
            if(l < dptr->ndet) continue; // a det in d didn't match grpdet[k]
/* d inherits the subgroup of grpdet[k] */
            dptr->sgrp = grpdet[k]->sgrp;
         }
      }
   }

   if(ngpd != NULL) free(ngpd);
   if(nalloc != NULL) free(nalloc);
   for(i=0; i<linkpar->ndet; i++) if(detgrp[i] != NULL) free(detgrp[i]);
   if(detgrp != NULL) free(detgrp);
   if(grpdet != NULL) free(grpdet);
   if(gsize != NULL) free(gsize);

#ifdef SUBGRPTEST
   printf("Done with group %d\n", qptr[0]->grp.det.dgrp);
   fflush(stdout);
#endif
   return;
}



/* Group an array of quad pointers using their fitted 6D posn */
/* (single friend, not FoF) */
int qgroupy(int n, QUAD **q, LINKPAR *linkpar, int *stack)
{
   int j, nstk, qidx, subgrp;
   double S[2*NDIM];
   void *kd, *set;

/* detection lists for each subgroup: depend on calloc to 0 everything */
   DETLIST *dg = (DETLIST *)calloc(n, sizeof(DETLIST));

/* Put all the unit vectors and angular velocities into a 6D kd tree */
   kd = kd_create(2*NDIM);
   
/* Insert all the unit vectors from these quads, with quad pointers */
   for(j=0; j<n; j++) {
      S[0] = q[j]->grp.X[0] / linkpar->subgrpxtol;
      S[1] = q[j]->grp.X[1] / linkpar->subgrpxtol;
      S[2] = q[j]->grp.X[2] / linkpar->subgrpxtol;
      S[3] = q[j]->grp.W[0] / linkpar->subgrpwtol;
      S[4] = q[j]->grp.W[1] / linkpar->subgrpwtol;
      S[5] = q[j]->grp.W[2] / linkpar->subgrpwtol;
      kd_insert(kd, S, (void *)(q+j));
   }

/* Loop over all points, assigning groups */
   subgrp = 0;
   for(j=0; j<n; j++) {         // j = index of array of quad pointers
/* Already grouped?  Continue until something not already grouped. */
      if(q[j]->grp.det.sgrp >= 0) continue;

/* New subgroup */
      subgrp += 1;
      nstk = 0;
      stack[nstk++] = j;
      q[j]->grp.det.sgrp = subgrp;
/* Number of detections encompassed by this subgroup's quads */
      dg[subgrp-1].ndet = 0;
      dg[subgrp-1].sgrp = subgrp;
      consistent_det(q[j], dg+subgrp-1);

/* What are the points within tol of the quad at the top of the stack? */
      S[0] = q[j]->grp.X[0] / linkpar->subgrpxtol;
      S[1] = q[j]->grp.X[1] / linkpar->subgrpxtol;
      S[2] = q[j]->grp.X[2] / linkpar->subgrpxtol;
      S[3] = q[j]->grp.W[0] / linkpar->subgrpwtol;
      S[4] = q[j]->grp.W[1] / linkpar->subgrpwtol;
      S[5] = q[j]->grp.W[2] / linkpar->subgrpwtol;
      set = kd_nearest_range(kd, S, 1.0);

//      printf("Nset= %d for quad %d\n", kd_res_size(set), j);

/* Examine each one, skipping if quad only matches itself */
      if(kd_res_size(set) > 1) {
         while( !kd_res_end(set) ) {
/* get the index in the quad pointer array of the current result item */
            qidx = (QUAD **)kd_res_item(set, NULL) - q;
#if 0
            printf("dgrp= %3d quad= %2d nstk= %2d thisubgrp= %2d idx= %8d qsubgrp= %2d chi= %8.3f %8.4f\n", q[j]->grp.det.dgrp, j, nstk, subgrp, qidx, q[qidx]->grp.det.sgrp, q[qidx]->chi, q[qidx]->grp.chin);
#endif

/* assign it to this group, skipping self */
            if(qidx != j && consistent_det(q[qidx], dg+subgrp-1)) {
               q[qidx]->grp.det.sgrp = subgrp;
               nstk++;
            }

/* go to the next entry */
            kd_res_next(set);
         }
      }
/* Free this result set */
      kd_res_free(set);

/* If a singleton, change group ID back to zero and decrement group number */
      if(nstk == 1) {
	 subgrp -= 1;
	 q[j]->grp.det.sgrp = 0;
      }
   }

/* Consolidate subgroups */
//   struct timeval tv0, tv1, tv2;
//   char info[1024];
//   double dt0, dt1;
//   gettimeofday(&tv0, NULL);
//   tv1 = tv0;
   consolidate_subgrps(linkpar, subgrp, dg, n, q);
//   sprintf(info, "dgrp %d sgrp %d with %d quads", q[0]->grp.det.dgrp, subgrp, n);
//   TOD(info);

   for(j=0; j<n; j++) if(dg[j].d != NULL) free(dg[j].d);
   if(dg != NULL) free(dg);

   kd_free(kd);
   return(subgrp);
}


/* Split each detection group into FoF subgroups with overlapping quads */
int quadgroup(int nquad, QUAD *quad, LINKPAR *linkpar)
{
   int i, n, k, maxgrp, nqok;

/* List of all quads */
   QUAD **q = (QUAD **)calloc(nquad, sizeof(QUAD*));
/* Groups numbers of all quads */
   int *grp = (int *)calloc(nquad, sizeof(int));

/* scratch stack for qgroupy */
   int *stack = (int *)calloc(nquad, sizeof(int));

   maxgrp = 0;

/* Make a list of all quads with their detection groups */
   for(k=nqok=0; k<nquad; k++) {
/* Inhibit bad quads */
      if(quad[k].grp.det.dgrp > 0) {
         grp[nqok] = quad[k].grp.det.dgrp;
         q[nqok] = quad + k;
         nqok++;
      }
   }
/* Sort list of quad pointers by detection group */
   tsort_ip(nqok, grp, (void **)q);

/* Visit each detection group and split into quad subgroups */
   for(k=0; k<nqok; k+=n) {
/* Identify runs of quads in the same detection group */
      for(n=1; n<nqok-k && grp[k+n] == grp[k+n-1]; n++);
      if(n == 1) {
         q[k]->grp.det.sgrp = 0;
      } else {
#ifdef SUBGRPTEST
         printf("Group %d\n", k);
#endif
         i = qgroupy(n, q+k, linkpar, stack);
         maxgrp = MAX(i, maxgrp);
      }
   }

   free(stack);
   free(q);
   free(grp);

   return(maxgrp);
}



/* Group all detections that are linked by a quad */
int detgroup(int nquad, QUAD *quad)
{
   int i, j, k, n, ng, idx, nnew, ig, n1, n2;
   int ndet, ngrp, nmax, *g;
   size_t nalloc;
   DET *d;
   DETLIST *grp;

/* Maximum count of detections in a quad */
   for(k=nmax=ndet=0; k<nquad; k++) {
      n = quad[k].v1->ndet + quad[k].v2->ndet;
      if(n > nmax) nmax = n;
      ndet += n;
   }   
   g = (int *)calloc(nmax, sizeof(int));

/* Allocate an initial group array big enough for any possible set of groups */
   grp = (DETLIST *)calloc((ndet+3)/4, sizeof(DETLIST));
   ngrp = 0;

/* Use each quad's linkage to group detections */
   for(k=0; k<nquad; k++) {

/* What are the detections of this quad? */
      n1 = quad[k].v1->ndet;
      n2 = quad[k].v2->ndet;

#ifdef GRPDEBUG
      printf("quad %d %s %s %s %s\n", k, d[0]->id, d[1]->id, d[2]->id, d[3]->id);
#endif

/* g[0:ng-1]=group IDs of all groups present among the detections */
      for(i=ng=0; i<n1+n2; i++) {
         if(i<n1) d = quad[k].v1->d[i];
         else     d = quad[k].v2->d[i-n1];
         if((g[ng] = d->dgrp) > 0) ng++;
      }

#ifdef GRPDEBUG
      printf("%d groups %d %d %d %d\n", ng, d[0]->dgrp, d[1]->dgrp, d[2]->dgrp, d[3]->dgrp);
#endif

/* Sort and remove dups: g[0:ng-1]=group IDs of unique groups, sorted */
      if(ng > 1) {
         tsort_i(ng, g, NULL);
         for(i=j=1; i<ng; i++) if(g[i] != g[j-1]) g[j++] = g[i];
         ng = j;
      }

#ifdef GRPDEBUG
      printf("%d unique groups %d %d %d %d\n", ng, g[0], g[1], g[2], g[3]);
#endif


/* No detection ever grouped before?  Create a new group. */
      if(ng == 0) {
/* Tell the detections which group they belong to.  NB: group n at index n-1 */
         nalloc = grp[ngrp].nalloc = grp[ngrp].ndet = n1 + n2;
         grp[ngrp].d = (DET **)calloc(nalloc, sizeof(DET *));
         for(i=0; i<nalloc; i++) {
            if(i<n1) d = quad[k].v1->d[i];
            else     d = quad[k].v2->d[i-n1];
            d->dgrp = ngrp+1;
            grp[ngrp].d[i] = d;
         }
         ngrp++;

#ifdef GRPDEBUG
         printf("new group %d\n", ngrp);
#endif

/* All detections that are grouped share the same group */
      } else if(ng == 1) {
/* Allocate space if required */
         idx = g[0] - 1;
         n = grp[idx].ndet;
         if(grp[idx].nalloc <= n + n1+n2 - 1) {
            grp[idx].nalloc += MAX(GRPALLOC, n1+n2);
            grp[idx].d = (DET **)realloc(grp[idx].d, 
                                    grp[idx].nalloc*sizeof(DET *));
#ifdef GRPDEBUG
            printf("realloc to %d n= %d\n", grp[idx].nalloc,grp[idx].ndet);
#endif
         }
/* Add unassigned detections to this group list, give them this group ID */
         for(i=0; i<n1+n2; i++) {
            if(i<n1) d = quad[k].v1->d[i];
            else     d = quad[k].v2->d[i-n1];
            if(d->dgrp < 0) {
               grp[idx].d[n++] = d;
               grp[idx].ndet += 1;
               d->dgrp = g[0];
            }
         }

#ifdef GRPDEBUG
         printf("add to group %d\n", g[0]);
#endif

/* More than one group present: merge */
      } else {
/* How many spaces do we need in the first group? */
         idx = g[0] - 1;
         n = grp[idx].ndet;
         nnew = n + n1+n2 - 2;
         for(i=1; i<ng; i++) nnew += grp[g[i]-1].ndet;
         if(grp[idx].nalloc <= nnew) {
            grp[idx].nalloc += MAX(nnew, GRPALLOC);
            grp[idx].d = (DET **)realloc(grp[idx].d, 
                                    grp[idx].nalloc*sizeof(DET *));
#ifdef GRPDEBUG
            printf("realloc to %d n= %d\n", grp[idx].nalloc,grp[idx].ndet);
#endif
         }

/* Add unassigned detections to the first group, give them this group ID */
         for(i=0; i<n1+n2; i++) {
            if(i<n1) d = quad[k].v1->d[i];
            else     d = quad[k].v2->d[i-n1];
            if(d->dgrp < 0) {
               grp[idx].d[n++] = d;
               grp[idx].ndet += 1;
               d->dgrp = g[0];
            }
         }

/* Merge the other groups onto the first, zero'ing their length */
         for(i=0; i<n1+n2; i++) {
            if(i<n1) d = quad[k].v1->d[i];
            else     d = quad[k].v2->d[i-n1];
            if( (ig = d->dgrp) == g[0]) continue;   // first group
            if(ig == -1) continue;                      // no group
            if(grp[ig-1].ndet == 0) continue;           // already wiped
#ifdef GRPDEBUG
            printf("Copy %d det from group %d to group %d with %d\n", grp[ig-1].ndet, ig, g[0], grp[idx].ndet);
#endif
            for(j=0; j<grp[ig-1].ndet; j++) {
               grp[idx].d[n++] = grp[ig-1].d[j];
               grp[idx].ndet += 1;
            }
            free(grp[ig-1].d);
            grp[ig-1].ndet = grp[ig-1].nalloc = 0;
#ifdef GRPDEBUG
            printf("WIPE OUT group %d\n", ig);
#endif
         }

/* Change the group IDs of all detections to the first */
         for(j=0; j<grp[idx].ndet; j++) grp[idx].d[j]->dgrp = g[0];

#ifdef GRPDEBUG
         printf("merge groups\n");
#endif
      }
   }

/* Squeeze out zero'ed groups: i=ok group, j=next slot */
   for(i=j=0; i<ngrp; i++) {
      if(grp[i].ndet == 0) continue;
      if(i > j) {
/* put group i in slot j */
#ifdef GRPDEBUG
  printf("Put group %d into slot %d\n", i+1, j+1);
#endif
         grp[j].ndet = grp[i].ndet;
         grp[j].nalloc = grp[i].nalloc;
         grp[j].d = grp[i].d;
/* Update group number for all detections in this group */
         for(k=0; k<grp[j].ndet; k++) grp[j].d[k]->dgrp = j+1;
/* Wipe out the group copy at slot i */
         free(grp[i].d);
         grp[i].ndet = grp[j].nalloc = 0;
      }
      j++;
   }
#ifdef GRPDEBUG
   printf("Squeeze from %d to %d\n", ngrp, j);
#endif
   ngrp = j;

/* Tell each quad what group it belongs to */
   for(k=0; k<nquad; k++) quad[k].grp.det.dgrp = quad[k].v1->d[0]->dgrp;

   free(g);
   free(grp);

   return(ngrp);
}


/* Return a sorted list of the unique detections in vectors v1 and v2 */
int unique_det(VEC *v1, VEC *v2, DET **d)
{
   int i, ndet, n1 = v1->ndet, n2 = v2->ndet;
   char **id = (char **)calloc(v1->ndet+v2->ndet, sizeof(char *));
   int *idx = (int *)calloc(v1->ndet+v2->ndet, sizeof(int));

   for(i=0; i<n1+n2; i++) {
      if(i<n1) d[i] = v1->d[i];
      else     d[i] = v2->d[i-n1];
      idx[i] = d[i]->idx;
   }
/* Sort on detection index (MJD,ID) */
   tsort_ip(n1+n2, idx, (void **)d);

/* Strip any dups */
   for(i=ndet=1; i<n1+n2; i++) if(idx[i] != idx[i-1]) d[ndet++] = d[i];

   free(id);
   free(idx);
   return(ndet);
}



#define QUADALLOC 1000

/* Match all vectors according to tol, test pumalink and puma quad */
int pumaquad(void *puma, int nvec, VEC *vec, LINKPAR *linkpar, QUAD **quad)
{
   int i, j, n, nquad, nalloc, detalloc, ndet;
   int dxrock=0, nchin=0, maxrock=0, maxmaxrock=0, nxrock=0;
   void *kd, *set;
   float pos1[2*NDIM], pos2[2*NDIM];
   double chi, chin, chin1[3], sw[5], rv[5], X[3], Wtan[3], dr=atan(1)/45;
   VEC *v1ptr;
   DET **d;

/* Detection lists for each quad */
   detalloc = QUADALLOC;
   d = (DET **)calloc(detalloc, sizeof(DET *));

   *quad = NULL;
   nquad = nalloc = 0;

   kd = kd_create(2*NDIM);
   
/* Insert all the entries from the first set of vectors into a kd-tree */
   for(j=0; j<nvec; j++) {
      if(vec[j].sigma <= 0) continue;          // skip failed eph
      if(linkpar->split && vec[j].d[0]->mjd > linkpar->mjd) continue;  // skip 2nd epoch
      for(i=0; i<2*NDIM; i++) pos1[i] = vec[j].S[i] / 
           ((2*(i/2) == i) ? linkpar->dxmax : linkpar->dwmax );
      kd_insertf(kd, pos1, (void *)&vec[j]);
   }

/* Match all the entries from the second set of vectors to the kd-tree */
   for(j=n=0; j<nvec; j++) {
      if(vec[j].sigma <= 0) continue;          // skip failed eph
/* Skip vectors which are in the first epoch */
      if(linkpar->split && vec[j].d[1]->mjd < linkpar->mjd) continue;

      for(i=0; i<2*NDIM; i++) pos2[i] = vec[j].S[i] / 
           ((2*(i/2) == i) ? linkpar->dxmax : linkpar->dwmax );

      set = kd_nearest_rangef(kd, pos2, 1.0);

      while( !kd_res_end(set) ) {

/* get the data and position of the current result item */
         v1ptr = (VEC *)kd_res_item(set, NULL);

/* Skip vectors which are not time ordered */
         if(vec[j].d[0]->mjd <= v1ptr->d[1]->mjd) {
            kd_res_next(set);   // go to the next entry
            continue;
         }

         if(TRACK_SINGLE_QUAD != NULL) {
            fprintf(TRACK_SINGLE_QUAD, "kdmch: %s %s %s %s %8d %8d\n",
                    v1ptr->d[0]->id, v1ptr->d[1]->id, 
                    vec[j].d[0]->id, vec[j].d[1]->id, j, n);
            fprintf(TRACK_SINGLE_QUAD, "S1: %8.4f %8.4f %8.4f %10.3e %10.3e %10.3e\n",
                    v1ptr->S[0], v1ptr->S[2], v1ptr->S[4],
                    v1ptr->S[1], v1ptr->S[3], v1ptr->S[5]);
            fprintf(TRACK_SINGLE_QUAD, "S2: %8.4f %8.4f %8.4f %10.3e %10.3e %10.3e\n",
                    vec[j].S[0], vec[j].S[2], vec[j].S[4],
                    vec[j].S[1], vec[j].S[3], vec[j].S[5]);
         }

/* Now do pumalink test, then pumaquad test, then keep if success */
         chi = pumalink(v1ptr, vec+j, linkpar->mjd, sw);
#ifdef LINK_DEBUG
         printf("chi= %10.3e\n", chi);
#endif
         if(TRACK_SINGLE_QUAD != NULL) {
            fprintf(TRACK_SINGLE_QUAD, "chi= %10.3e\n", chi);
         }

         if(chi <= linkpar->chimax && chi > 0) {
            nxrock = NXROCK;

            if(v1ptr->ndet + vec[j].ndet > detalloc) {
               detalloc += v1ptr->ndet + vec[j].ndet;
               d = (DET **)realloc(d, detalloc*sizeof(DET *));
            }
/* Get the unique detections from these two vectors */
            ndet = unique_det(v1ptr, vec+j, d);
/* Do a puma fit to all the detections */
            chin = pumatest(puma, ndet, d, linkpar, X, Wtan, rv, chin1);
#ifdef LINK_DEBUG
            printf("chin= %10.3e\n", chin);
#endif
            if(TRACK_SINGLE_QUAD != NULL) {
               fprintf(TRACK_SINGLE_QUAD, "chin= %10.3e\n", chin);
            }

            if(chin <= linkpar->chinmax && chin > 0 && chin*chi < linkpar->chiprod) {
/* Save this quad */
               if(nquad >= nalloc) {
                  *quad = (QUAD *)realloc(*quad, (nalloc+QUADALLOC)*sizeof(QUAD));
                  nalloc += QUADALLOC;
               }
#ifdef LINK_DEBUG
               printf("nquad= %d\n", nquad);
#endif
               if(TRACK_SINGLE_QUAD != NULL) {
                  fprintf(TRACK_SINGLE_QUAD, "nquad= %d\n", nquad);
               }
/* Save the pumalink parameters as part of this quad */
               (*quad)[nquad].v1 = v1ptr;
               (*quad)[nquad].v2 = vec + j;
               (*quad)[nquad].chi = LOGCLIP(chi);
               (*quad)[nquad].s0 = sw[0];
               (*quad)[nquad].w0 = sw[1];
               (*quad)[nquad].ds = sw[2];
               (*quad)[nquad].dw = sw[3];
               (*quad)[nquad].swcov = sw[4];
/* Save the puma fit parameters as this quad's grp */
               (*quad)[nquad].grp.chin = LOGCLIP(chin);
               (*quad)[nquad].grp.chin1 = LOGCLIP(chin1[0]);
               (*quad)[nquad].grp.chin2 = LOGCLIP(chin1[1]);
               (*quad)[nquad].grp.rfit = rv[0];
               (*quad)[nquad].grp.wfit = rv[1] / rv[0];
               (*quad)[nquad].grp.dlnr = rv[2];
               (*quad)[nquad].grp.dwr = rv[3];
               (*quad)[nquad].grp.rwcov = rv[4];
               (*quad)[nquad].grp.X[0] = X[0];
               (*quad)[nquad].grp.X[1] = X[1];
               (*quad)[nquad].grp.X[2] = X[2];
               (*quad)[nquad].grp.W[0] = Wtan[0] / dr * SECDAY;
               (*quad)[nquad].grp.W[1] = Wtan[1] / dr * SECDAY;
               (*quad)[nquad].grp.W[2] = Wtan[2] / dr * SECDAY;
               (*quad)[nquad].grp.det.dgrp = -1;
               (*quad)[nquad].grp.det.sgrp = -1;
/* Save all the detections as part of the grp */
               (*quad)[nquad].grp.det.ndet = ndet;
               (*quad)[nquad].grp.det.nalloc = ndet;
               (*quad)[nquad].grp.det.d = (DET **)calloc(ndet, sizeof(DET *));
               for(i=0; i<ndet; i++) (*quad)[nquad].grp.det.d[i] = d[i];

               if(DEBUG_QUAD != NULL) {
                  fprintf(DEBUG_QUAD, "%s %s %s %s %4d.%-3d %6.3f %6.3f %6.3f %6.3f %7.4f %8.5f %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f\n",
                          (*quad)[nquad].v1->d[0]->id, (*quad)[nquad].v1->d[1]->id, 
                          (*quad)[nquad].v2->d[0]->id, (*quad)[nquad].v2->d[1]->id, 
                          (*quad)[nquad].grp.det.dgrp, (*quad)[nquad].grp.det.sgrp,
                          (*quad)[nquad].chi, (*quad)[nquad].grp.chin,
                          (*quad)[nquad].grp.chin1, (*quad)[nquad].grp.chin2,
                          (*quad)[nquad].grp.rfit/AU, (*quad)[nquad].grp.wfit*SECDAY,
                          (*quad)[nquad].grp.X[0], (*quad)[nquad].grp.X[1], (*quad)[nquad].grp.X[2],
                          (*quad)[nquad].grp.W[0], (*quad)[nquad].grp.W[1], (*quad)[nquad].grp.W[2]);
               }

#if 0
               printf("%s %s %s %s  %8.2f  %9.4f  %9.4f %9.4f  %9.4f %9.4f %9.4f %9.4f %10.3e %10.3e %10.3e\n",
                      v1ptr->d[0]->id, v1ptr->d[1]->id, 
                      vec[j].d[0]->id, vec[j].d[1]->id, chi, chin,
                      chin1[0], chin1[1],
                      rv[0]/AU, rv[1]/rv[0]*SECDAY,
                      X[0], X[1], X[2],
                      sqrt(Wtan[0]*Wtan[0]+Wtan[1]*Wtan[1]+Wtan[2]*Wtan[2])/dr*SECDAY, Wtan[0], Wtan[1], Wtan[2]);
#endif
               nquad++;
            }
            maxrock = MAX(maxrock, NXROCK-nxrock);
            maxmaxrock = MAX(maxmaxrock, NXROCK-nxrock);
            nchin++;
         }

/* go to the next entry */
         kd_res_next(set);
         n++;
         if((n%10000)==0) {
            fprintf(stderr, "Quad %8d %8d %8d %8d %8.1f %8d %12d   \r",
                    n, nquad, maxrock, maxmaxrock,
                    ((double)(NXROCK-dxrock))/MAX(1,nchin), NXROCK-dxrock, 
                    NXROCK);
            nchin = 0;
            maxrock = 0;
            dxrock = NXROCK;
         }
         
      }
      kd_res_free(set);
   }
   fprintf(stderr, "\n");

   kd_free(kd);
   free(d);

   return(nquad);
}

/* Return chin from a puma fit to all detections of a pair of vectors */
double pumatest(void *puma, int ndet, DET **d, LINKPAR *linkpar,
                double *X, double *Wtan, double *rv, double *chin1)
{
   int i;
   double mjd, sigma, chin, lnr, wr, dlnr, dwr, rwcov, ephmjd;
   double V[3], dX[3], dV[3], rho, rhodot;

   ephmjd = linkpar->mjd;

   puma_reset(puma);

/* Load all the detections into puma */
   for(i=0; i<ndet; i++) {
      puma_obs(puma, 0, d[i]->mjd, 
               d[i]->ra, d[i]->dec, (double)d[i]->xerr, (double)d[i]->terr,
               (double)d[i]->lng, (double)d[i]->lat, (double)d[i]->elev);
   }

/* Fit for lnr, wr */
   chin = puma_gridfit(puma, &lnr, &wr, &dlnr, &dwr, &rwcov, linkpar);
   chin1[0] = puma_const(puma, chin1+1);
//   double rfit, vfit;
//   rfit = exp(lnr);                             // [AU]
//   vfit = wr / SECDAY * rfit*AU * 1e-3;         // [km/s]

/* Evaluate barycenter location at ephmjd */
   i = puma_obs(puma, 1, ephmjd, 0., 0., 0., 0., 0., 0., 0.);
   puma_bary(puma, i, linkpar->integrate, &mjd, X, V, dX, dV, &sigma);

/* At the requested epoch... */
   rho = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
   rhodot = (X[0]*V[0]+X[1]*V[1]+X[2]*V[2]) / rho;

/* Return particulars of the fit */
   rv[0] = rho;      // [m]
   rv[1] = rhodot;   // [m/s]
   rv[2] = dlnr;
   rv[3] = dwr;
   rv[4] = rwcov;

/* Return unit vector */
   X[0] /= rho;
   X[1] /= rho;
   X[2] /= rho;

/* Return Vtan = tangential angular velocity [rad/sec] */
   Wtan[0] = (V[0] - X[0] * rhodot) / rho;
   Wtan[1] = (V[1] - X[1] * rhodot) / rho;
   Wtan[2] = (V[2] - X[2] * rhodot) / rho;

   return(chin);
}


/* InvCov product of two vectors: v1^T C v2;  C=(c0 c1 \ c1 c2) */
#define CDOT(v1,C,v2) (v1[0]*(C[0]*v2[0]+C[1]*v2[1])+v1[1]*(C[1]*v2[0]+C[2]*v2[1]))

/* InvCov product of two 2D vectors: (v1-z1)^T C (v2-z2); C=(c0 c1 \ c1 c2) */
#define XDOT(v1,z1,C,v2,z2) ( (*(v1)-*(z1)) * (C[0] * (*(v2)-*(z2)) + C[1] * (*(v2+1)-*(z2+1))) + (*(v1+1)-*(z1+1)) * (C[1] * (*(v2)-*(z2)) + C[2]*(*(v2+1)-*(z2+1))) )

/* Return chi^2 from a pumalink assessment of a pair of vectors */
double pumalink(VEC *v1, VEC *v2, double ephmjd, double *sw)
{
   int i, j, k;

   double dt1, Dt1, dt2, Dt2, det, A, B;;

/* Covariance matrices: [3] for 2x2 symmetric matrix A11,A12,A22 */
   double C1[NDIM][3], C2[NDIM][3];
   double C1inv[NDIM][3], C2inv[NDIM][3];
   double C[NDIM][3];
/* Polynomial residuals */
   double X1[NPAR][2*NDIM], X2[NPAR][2*NDIM], X0[NPAR][2*NDIM];
#ifdef BUDGE_CHI_BY_DSW
   double dX0[NPAR][2*NDIM], dchi[NPAR*NPAR];  /* (1),ds,dw linear terms */
#endif

/* chi^2 polynomial terms from squared fits */
   double chi[NPAR*NPAR];

   dt1 = v1->dt;
   dt2 = v2->dt;

   Dt1 = (ephmjd - v1->mjd) * SECDAY;
   Dt2 = (ephmjd - v2->mjd) * SECDAY;

/* Reorganize float fit coefficients into double arrays */
   for(k=0; k<2*NDIM; k++) {
      for(i=0; i<NPAR; i++) {
         X1[i][k] = v1->A[i+k*NPAR];
         X2[i][k] = v2->A[i+k*NPAR];
      }
   }

#ifdef LINK_DEBUG
   for(k=0; k<2*NDIM; k++) printf("vecs: %10.3e %10.3e %10.3e   %10.3e %10.3e %10.3e   %10.3e   %10.3e %10.3e\n", X1[0][k], X1[1][k], X1[2][k], X2[0][k], X2[1][k], X2[2][k], X1[0][k]-X2[0][k], v1->sigma, v2->sigma);
#endif

/* Calculate covariance matrices */
   for(k=0; k<NDIM; k++) {
/* First 2x2 vector covariance for x,vx from point uncertainty */
      C1[k][0] = 2 * VARBLOAT * (v1->sigma * v1->sigma) / (dt1*dt1) * (Dt1*Dt1+dt1*dt1/4);
      C1[k][1] = 2 * VARBLOAT * (v1->sigma * v1->sigma) / (dt1*dt1) * Dt1;
      C1[k][2] = 2 * VARBLOAT * (v1->sigma * v1->sigma) / (dt1*dt1);

#ifdef LINK_DEBUG
//      double ctest=C1[k][1]/sqrt(C1[k][0]*C1[k][2]);
//      if(ctest >  MAXCRV) ctest =  MAXCRV;
//      if(ctest < -MAXCRV) ctest = -MAXCRV;
//      C1[k][1] = ctest * sqrt(C1[k][0]*C1[k][2]);

      printf("rms: %10.3e %10.3e %6.3f", v1->rms[0+2*k], v1->rms[1+2*k], v1->cov[k]);
      printf("  cov1: %10.3e %6.3f %10.3e", sqrt(C1[k][0]), C1[k][1]/sqrt(C1[k][0]*C1[k][2]), sqrt(C1[k][2]));
#endif

/* Add s,w fit covariance */
      if(v1->cov[k] >  MAXCOV) v1->cov[k] =  MAXCOV; // sanity
      if(v1->cov[k] < -MAXCOV) v1->cov[k] = -MAXCOV;
      C1[k][0] += v1->rms[0+2*k] * v1->rms[0+2*k];
      C1[k][2] += v1->rms[1+2*k] * v1->rms[1+2*k];
      C1[k][1] += v1->rms[0+2*k] * v1->rms[1+2*k] * v1->cov[k];

#ifdef LINK_DEBUG
      printf("  %10.3e %6.3f %10.3e", sqrt(C1[k][0]), C1[k][1]/sqrt(C1[k][0]*C1[k][2]), sqrt(C1[k][2]));
#endif

/* Invert */
      det = C1[k][0]*C1[k][2] - C1[k][1]*C1[k][1];
      if(det <= 0) {     // give up on s,w covariance, should never occur
         C1[k][1] += v1->rms[0+2*k] * v1->rms[1+2*k] * v1->cov[k];
         det = C1[k][0]*C1[k][2] - C1[k][1]*C1[k][1];
      }
      C1inv[k][0] =  C1[k][2] / det;
      C1inv[k][1] = -C1[k][1] / det;
      C1inv[k][2] =  C1[k][0] / det;

#ifdef LINK_DEBUG
      printf("  %10.3e %6.3f %10.3e\n", sqrt(C1inv[k][0]), C1inv[k][1]/sqrt(C1inv[k][0]*C1inv[k][2]), sqrt(C1inv[k][2]));
#endif

/* Second vector covariance for x,vx from point uncertainty */
      C2[k][0] = 2 * VARBLOAT * (v2->sigma * v2->sigma) / (dt2*dt2) * (Dt2*Dt2+dt2*dt2/4);
      C2[k][1] = 2 * VARBLOAT * (v2->sigma * v2->sigma) / (dt2*dt2) * Dt2;
      C2[k][2] = 2 * VARBLOAT * (v2->sigma * v2->sigma) / (dt2*dt2);

#ifdef LINK_DEBUG
//      ctest = C2[k][1]/sqrt(C2[k][0]*C2[k][2]);
//      if(ctest >  MAXCRV) ctest =  MAXCRV;
//      if(ctest < -MAXCRV) ctest = -MAXCRV;
//      C2[k][1] = ctest * sqrt(C2[k][0]*C2[k][2]);
      printf("rms: %10.3e %10.3e %6.3f", v2->rms[0+2*k], v2->rms[1+2*k], v2->cov[k]);
      printf("  cov2: %10.3e %6.3f %10.3e", sqrt(C2[k][0]), C2[k][1]/sqrt(C2[k][0]*C2[k][2]), sqrt(C2[k][2]));
#endif

/* Add s,w fit covariance */
      if(v2->cov[k] >  MAXCOV) v2->cov[k] =  MAXCOV; // sanity
      if(v2->cov[k] < -MAXCOV) v2->cov[k] = -MAXCOV;
      C2[k][0] += v2->rms[0+2*k] * v2->rms[0+2*k];
      C2[k][2] += v2->rms[1+2*k] * v2->rms[1+2*k];
      C2[k][1] += v2->rms[0+2*k] * v2->rms[1+2*k] * v2->cov[k];

#ifdef LINK_DEBUG
      printf("  %10.3e %6.3f %10.3e", sqrt(C2[k][0]), C2[k][1]/sqrt(C2[k][0]*C2[k][2]), sqrt(C1[k][2]));
#endif

/* Invert */
      det = C2[k][0]*C2[k][2] - C2[k][1]*C2[k][1];
      if(det <= 0) {     // give up on s,w covariance, should never occur
         C2[k][1] -= v2->rms[0+2*k] * v2->rms[1+2*k] * v2->cov[k];
         det = C2[k][0]*C2[k][2] - C2[k][1]*C2[k][1];
      }
      C2inv[k][0] =  C2[k][2] / det;
      C2inv[k][1] = -C2[k][1] / det;
      C2inv[k][2] =  C2[k][0] / det;

#ifdef LINK_DEBUG
      printf("  %10.3e %6.3f %10.3e\n", sqrt(C2inv[k][0]), C2inv[k][1]/sqrt(C2inv[k][0]*C2inv[k][2]), sqrt(C2inv[k][2]));
#endif

/* Inverse of sum of inverse covariance for each coordinate */
      det=(C1inv[k][0] + C2inv[k][0]) * (C1inv[k][2] + C2inv[k][2]) -
          (C1inv[k][1] + C2inv[k][1]) * (C1inv[k][1] + C2inv[k][1]);
      C[k][0] =  (C1inv[k][2] + C2inv[k][2]) / det;
      C[k][1] = -(C1inv[k][1] + C2inv[k][1]) / det;
      C[k][2] =  (C1inv[k][0] + C2inv[k][0]) / det;
   }



/* Evaluate Smin(s,w) = C(Cinv1 S1 + Cinv2 S2) that minimizes chi^2 */
   for(k=0; k<NDIM; k++) {
/* Sum of coefficients for each fit term */
      for(j=0; j<NPAR; j++) {
         A = C1inv[k][0] * X1[j][0+2*k] + C1inv[k][1] * X1[j][1+2*k];
         B = C1inv[k][1] * X1[j][0+2*k] + C1inv[k][2] * X1[j][1+2*k];
#ifdef BUDGE_CHI_BY_DSW
/* ds, dw terms on the S1 variable */
         dX0[j][0+2*k] = C[k][0]*A + C[k][1]*B;
         dX0[j][1+2*k] = C[k][1]*A + C[k][2]*B;
#endif
/* Full 1,s,w,... coefficients for Smin */
         A += C2inv[k][0] * X2[j][0+2*k] + C2inv[k][1] * X2[j][1+2*k];
         B += C2inv[k][1] * X2[j][0+2*k] + C2inv[k][2] * X2[j][1+2*k];
         X0[j][0+2*k] = C[k][0]*A + C[k][1]*B;
         X0[j][1+2*k] = C[k][1]*A + C[k][2]*B;
      }
   }
// Smin(1,s,w,... ds,dw) for each state vector component
// X0[j][k] = coeffs of Smin[j:1,s,w...] for [k:x,vx, y,vy, z,vz]
// dX0[j][k] = coeffs of Smin[j:-,ds,dw...] for [k:x,vx, y,vy, z,vz]


#ifdef LINK_DEBUG
   for(k=0; k<2*NDIM; k++) printf("Smin: %10.3e  %10.3e  %10.3e  %10.3e  %10.3e\n", X0[0][0+k], X0[1][0+k], X0[2][0+k], X1[0][k]-X0[0][k], X2[0][k]-X0[0][k]);
#endif


/* chi^2 1,s,w... coefficients: save upper diagonal symmetric matrix */
   for(i=0; i<NPAR*NPAR; i++) chi[i] = 0;
#ifdef BUDGE_CHI_BY_DSW
   for(i=0; i<NPAR*NPAR; i++) dchi[i] = 0;
#endif

// chi^2(s,w) = (S1-Smin)T C1inv (S1-Smin) + (S2-Smin)T C2inv (S2-Smin)
// The covariance matrix is 3 blocks of 2x2 for x,vx  y,vy  z,vz
// Therefore do the matrix multiplication 
// Coefficients laid out in a partially filled matrix as
//                  1  s    w
//  chi[i+j*NPAR] = .  s^2  sw
//                  .  .    w^2

/* Add the coefficients for the 6 chi^2 terms as a symmetric matrix: 1 s w */
   for(k=0; k<NDIM; k++) {
      for(j=0; j<NPAR; j++) {    // chi[ 0   1   2   3   4   5   6   7   8 ]
         for(i=j; i<NPAR; i++) { //      1   s   w   .  s^2  sw  .   .  w^2
            chi[i+j*NPAR] += ((i==j)?1:2)*XDOT(X1[i]+2*k, X0[i]+2*k, C1inv[k],
                                               X1[j]+2*k, X0[j]+2*k);
            chi[i+j*NPAR] += ((i==j)?1:2)*XDOT(X2[i]+2*k, X0[i]+2*k, C2inv[k],
                                               X2[j]+2*k, X0[j]+2*k);
         }
      }
   }
   
/* Solve for chi^2 = chi0 + e1(s-s0^2) + e2(s-s0)(w-w0) + e3(w-w0)^2 */
   double s0, w0, chi0, ds, dw, csw;
   det = 4*chi[4]*chi[8] - chi[5]*chi[5];
   s0 = (chi[5]*chi[2] - 2*chi[8]*chi[1]) / det;
   w0 = (chi[5]*chi[1] - 2*chi[4]*chi[2]) / det;
   chi0 = chi[0] - chi[4]*s0*s0 - chi[5]*s0*w0 - chi[8]*w0*w0;

/* Invert the chi^2 quadratic terms to the covariance matrix */
   ds = 4*chi[8] / det;
   dw = 4*chi[4] / det;

   csw = -4*chi[5] / det;
   if(ds > 0) ds = sqrt(ds);
   else       ds = -sqrt(-ds);
   if(dw > 0) dw = sqrt(dw);
   else       dw = -sqrt(-dw);
   csw /= ds * dw;

   sw[0] = s0;
   sw[1] = w0;
   sw[2] = ds;
   sw[3] = dw;
   sw[4] = csw;
   return(chi0);
}


/* Project one vector to the ephemeris time */
int pumephem(LINKPAR *linkpar, void *puma, VEC *vec, FILE *fpscat, int nscat)
{
   int i, j, k, l, m, npt, nvgrid, ephdet;

   double R, VR, Dt, chin, chinsum, lnr, wr, mjd, RMS[4], v2;
   double smax, smin, stest, vmin, vmax;
   double dr=atan(1)/45;
   double dX[NDIM], dV[NDIM], S[2*NDIM];
   double X[MAXPAR+7/*plenty*/], V[MAXPAR+7/*plenty*/];
   double DTMP[NPAR*6], DX[NPAR*6], DD[2*NDIM], DV[NDIM];
   double A[NPAR*NDIM*2];
   double M[NPAR*NPAR], MINV[NPAR*NPAR];

   double s, w, rho, rhodot, sigma, vtan, ang;
//   char varname[MAXPAR][10], coordname[2*NDIM][10];

/* Signal failed fit by setting vec->sigma = -1 */
   vec->sigma = -1;

/* Initialize linear fit */
   npt = 0;
   for(k=0; k<2*NDIM; k++) DD[k] = 0;
   for(k=0; k<NDIM; k++) DV[k] = 0;
   for(i=0; i<NPAR*6; i++) DX[i] = 0;
   for(i=0; i<NPAR*NPAR; i++) M[i] = 0;

/* Time interval from mid-vector to ephemeris [sec] */
   Dt = (linkpar->mjd - vec->mjd) * SECDAY;

   if(TRACK_SINGLE_VEC != NULL) {
      fprintf(TRACK_SINGLE_VEC, "enter: %s %s %7.1f\n",
              vec->d[0]->id, vec->d[1]->id, (vec->d[1]->mjd-vec->d[0]->mjd)*SECDAY);
   }

/* Give puma this vector's detections and the ephemeris time */
   puma_reset(puma);

   for(i=0; i<vec->ndet; i++) {
      puma_obs(puma, 0, vec->d[i]->mjd, vec->d[i]->ra, vec->d[i]->dec,
               (double)vec->d[i]->xerr, (double)vec->d[i]->terr,
               (double)vec->d[i]->lng, (double)vec->d[i]->lat,
               (double)vec->d[i]->elev);
//      printf("%10.4f %8.4f %8.4f %s\n", vec->d[i]->mjd,vec->d[i]->ra, vec->d[i]->dec, vec->d[i]->id);
   }

/* Give puma the reference time and save the returned index */
   ephdet = puma_obs(puma, 1, linkpar->mjd, 0., 0., 0., 0., 0., 0., 0.);

/* Grid ranges */
   nvgrid = linkpar->nv;
   rho = rhodot = 0;

/* Fewer than 4 detections means use the default r,vr grid */
   if(vec->ndet < NEPHDETMIN) {
/* Step over the radius and velocity grid */
      smax = 1 / linkpar->rmin;
      smin = 1 / linkpar->rmax;

/* Sanity check these ranges, given the angular velocity and the Dt */
      ang = sphdist(vec->d[0]->ra, vec->d[0]->dec,
                    vec->d[1]->ra, vec->d[1]->dec);

/* Angular velocity [rad/sec] */
      vtan = ang / ABS(vec->d[1]->mjd - vec->d[0]->mjd) * dr/SECDAY;

/* Maximum r = minimum s with an acceptable physical velocity */
      stest = (vtan / (linkpar->vtmax*1e3)) * AU;
      if(stest > 1/linkpar->rmax) smin = stest;
      if(smin > 1/linkpar->rmin) smin = 0.5 / linkpar->rmin;

/* Maximum inward radial velocity consistet with the wDt limit */
      vmax = linkpar->vmax;
      vmin = linkpar->vmin;
//   if(vmin < -AU / (smin*Dt)) vmin = -AU / (smin*Dt);

// FIXME: not well thought through, basically disabled for now
/* NEPHDETMIN or more detections requires the r,vr grid to be consistent */
   } else {
      double dlnr, dwr, rwcov, GRIDSIG=3;
      chin = puma_gridfit(puma, &lnr, &wr, &dlnr, &dwr, &rwcov, linkpar);
      smin = 1 / exp(lnr+GRIDSIG*dlnr);         // [1/AU]
      smax = 1 / exp(lnr-GRIDSIG*dlnr);
      vmin = (wr-GRIDSIG*dwr) * exp(lnr)*AU / (1e3*SECDAY);     // [km/s]
      vmax = (wr+GRIDSIG*dwr) * exp(lnr)*AU / (1e3*SECDAY);
      if(vmax > linkpar->vmax) vmax = linkpar->vmax;
      if(vmin < linkpar->vmin) vmin = linkpar->vmin;
//      printf("gridfit: %8.4f %8.4f %8.4f %8.4f %12.2f\n", lnr, dlnr, wr, dwr, chin);
   }
   if(TRACK_SINGLE_VEC != NULL) {
      fprintf(TRACK_SINGLE_VEC, "%8.4f %8.4f %8.1f %8.1f\n", 1/smax, 1/smin, vmin, vmax);
   }

   for(k=0; k<2*NDIM; k++) vec->S[k] = 0;
   chinsum = 0;

/* Step in s=1/r */
   for(i=0; i<linkpar->nr; i++) {

/* Step over the velocity grid */
      if(linkpar->nv == 0) nvgrid = i + 1;        // triangular grid
      for(j=0; j<nvgrid; j++) {

//      R = 1 / (smin + (i * (smax-smin)) / MAX(1,linkpar->nr-1));
/* Step in lnr [AU] */
         R = exp((i*log(smin/smax))/MAX(1,linkpar->nr-1)) / smin;

/* Step in VR [km/s]*/
         if(nvgrid == 1) {
            VR = 0.5 * (vmin + vmax);
         } else {
            VR = vmin + (j * (vmax-vmin)) / (nvgrid-1);
         }

/* Tweak R so that it will end up closer to the desired R at ephmjd */
         if(vec->ndet < NEPHDETMIN) {
            R -= VR*1e3 * Dt / AU;
            if(R < RMIN) R = RMIN;   // Floor for R
         }

         lnr = log(R);
         wr = VR*1e3/(R*AU)*SECDAY;

/* Where does this end up? */
         chin = puma_chi(puma, lnr, wr, RMS);

/* Does this pass the chin test? */
         if(chin > linkpar->chieph) {
#ifdef EXTRAPOLATE_TEST
            dV[0] = V[0] - X[0]/rho * rhodot;
            dV[1] = V[1] - X[1]/rho * rhodot;
            dV[2] = V[2] - X[2]/rho * rhodot;
            vtan = sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]);
            ang = sphdist(vec->d[0]->ra, vec->d[0]->dec,
                          vec->d[1]->ra, vec->d[1]->dec);
            printf("Clip w= %.2f Dt= %.1f npt= %d rho,rhodot= %10.6f %10.3f R,VR=  %10.6f %10.3f Vtan= %10.3f ang= %8.3f\n",
                    (rhodot/1e3)/(rho/AU), Dt, npt, rho/AU, rhodot/1e3, R, VR, vtan/1e3, ang);
#endif
            if(TRACK_SINGLE_VEC != NULL) {
               fprintf(TRACK_SINGLE_VEC, "chi: %s %s (%d) %8.4f %10.3e %9.1f\n",
                       vec->d[0]->id, vec->d[1]->id, vec->ndet, R, VR, chin);
            }
            continue;
         }

/* Get the ecliptic location with respect to the barycenter */
         puma_bary(puma, ephdet, linkpar->integrate, &mjd, X, V, dX, dV, &sigma);

/* Bary-obj distance */
         rho = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
         rhodot = (X[0]*V[0]+X[1]*V[1]+X[2]*V[2]) / rho;

/* Does this pass the wDt test?  |rhodot[m/s]/rho[m] * Dt[sec]| < wdtmin */
         if(ABS(rhodot/rho*Dt) > linkpar->wdtmin) {
#ifdef EXTRAPOLATE_TEST
            dV[0] = V[0] - X[0]/rho * rhodot;
            dV[1] = V[1] - X[1]/rho * rhodot;
            dV[2] = V[2] - X[2]/rho * rhodot;
            vtan = sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]);
            ang = sphdist(vec->d[0]->ra, vec->d[0]->dec,
                          vec->d[1]->ra, vec->d[1]->dec);
            printf("Clip w= %.2f Dt= %.1f npt= %d rho,rhodot= %10.6f %10.3f R,VR=  %10.6f %10.3f Vtan= %10.3f ang= %8.3f\n",
                    (rhodot/1e3)/(rho/AU), Dt, npt, rho/AU, rhodot/1e3, R, VR, vtan/1e3, ang);
#endif
            if(TRACK_SINGLE_VEC != NULL) {
               fprintf(TRACK_SINGLE_VEC, "wdt: %s %s %7.1f %12.5e\n",
                       vec->d[0]->id, vec->d[1]->id, (vec->d[1]->mjd-vec->d[0]->mjd)*SECDAY, sigma);
            }
            continue;
         }

/* Does this pass the physical vmax test?  |V|^2 < vtmax[km/s]^2 */
         if( (v2=(V[0]*V[0]+V[1]*V[1]+V[2]*V[2])) > linkpar->vtmax*linkpar->vtmax*1e6) {
#ifdef EXTRAPOLATE_TEST
            dV[0] = V[0] - X[0]/rho * rhodot;
            dV[1] = V[1] - X[1]/rho * rhodot;
            dV[2] = V[2] - X[2]/rho * rhodot;
            ang = sphdist(vec->d[0]->ra, vec->d[0]->dec,
                          vec->d[1]->ra, vec->d[1]->dec);
            vtan = sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]);
            printf("Clip: V= %.2f Dt= %.1f npt= %d rho,rhodot= %10.6f %10.3f R,VR=  %10.6f %10.3f Vtan= %10.3f ang= %8.3f %8.3f\n",
                    sqrt(v2)/1e3, Dt, npt, rho/AU, rhodot/1e3, R, VR, vtan/1e3, ang, ang/(vec->d[1]->mjd-vec->d[0]->mjd));

            puma_bary(puma, 1, linkpar->integrate, &mjd, X, V, dX, dV, &sigma);
            printf("Puma1: %11.5f %9.5f %9.5f %9.5f  %9.2f %9.2f %9.2f\n",
                    mjd, X[0]/AU, X[1]/AU, X[2]/AU, V[0]/1e3, V[1]/1e3, V[2]/1e3);
            puma_bary(puma, 2, linkpar->integrate, &mjd, X, V, dX, dV, &sigma);
            printf("Puma2: %11.5f %9.5f %9.5f %9.5f  %9.2f %9.2f %9.2f\n",
                    mjd, X[0]/AU, X[1]/AU, X[2]/AU, V[0]/1e3, V[1]/1e3, V[2]/1e3);
#endif
            if(TRACK_SINGLE_VEC != NULL) {
               fprintf(TRACK_SINGLE_VEC, "vmax: %s %s %7.1f %12.5e\n",
                       vec->d[0]->id, vec->d[1]->id, (vec->d[1]->mjd-vec->d[0]->mjd)*SECDAY, sigma);
            }
            continue;
         }

/* Interleave x,v for the 3 coordinates, divide by distance */
         S[0] = X[0] / rho;  // r=infinity unit vector
         S[1] = V[0] / rho;
         S[2] = X[1] / rho;
         S[3] = V[1] / rho;  // r=inf tang ang vel
         S[4] = X[2] / rho;
         S[5] = V[2] / rho;
/* Contribute to an average location */
         for(k=0; k<2*NDIM; k++) vec->S[k] += S[k] / MAX(0.1, chin);
         chinsum += 1 / MAX(0.1, chin);

//         printf("%8.4f %7.1f %9.2f %8.5f %8.5f %8.5f %10.3e %10.3e %10.3e\n",                R, VR, chin, S[0], S[2], S[4], S[1], S[3], S[5]);

/* Sum linear fit contributions */
         s = 1 / (rho/AU);             // [AU^-1]     1/r
         w = (rhodot/1e3) / (rho/AU);  // [km/s/AU]  vr/r

/* What are the linear variables?  (up to MAXPAR=6) */
/* Note: overload of X[] as position variable for puma and state vector here */
         X[0] = 1;
         X[1] = s;
         X[2] = w;
         X[3] = w*w;
         X[4] = s*w;
         X[5] = s*s;

/* Add their contributions to data vector and matrix */
         for(l=0; l<NPAR; l++) {
            for(m=0; m<NPAR; m++) M[l+m*NPAR] += X[l] * X[m];
         }
         for(k=0; k<2*NDIM; k++) {
            DD[k] += S[k] * S[k];                  // sq for RMS
            if((k%2)==1) DV[k/2] += S[k-1] * S[k]; // x,v covariance
            // weighted state vector:
            for(l=0; l<NPAR; l++) DX[l+k*NPAR] += X[l] * S[k];
         }
         npt++;
      } // nvgrid
   } // nrgrid

   if(npt < NPAR) {
#ifdef BLAB_FIT_FAIL
      dV[0] = V[0] - X[0]/rho * rhodot;
      dV[1] = V[1] - X[1]/rho * rhodot;
      dV[2] = V[2] - X[2]/rho * rhodot;
      vtan = sqrt(dV[0]*dV[0]+dV[1]*dV[1]+dV[2]*dV[2]);
      ang = sphdist(vec->d[0]->ra, vec->d[0]->dec,
                    vec->d[1]->ra, vec->d[1]->dec);
      fprintf(stderr, "npt= %d < npar= %d Dt= %.1f Vtan= %10.3f ang= %8.3f %10.2f rho,rhodot= %10.6f %10.3f\n", 
              npt, NPAR, Dt,
              vtan/1e3, ang, ang/(vec->d[0]->mjd-vec->d[1]->mjd),
              rho/AU, rhodot/1e3);
#endif
      if(TRACK_SINGLE_VEC != NULL) {
         fprintf(TRACK_SINGLE_VEC, "npt: %s %s %7.1f %12.5e\n",
                 vec->d[0]->id, vec->d[1]->id, (vec->d[1]->mjd-vec->d[0]->mjd)*SECDAY, sigma);
      }
      return(-2);
   }

/* Solve for linear fits, and write the results */
   memcpy(MINV, M, NPAR*NPAR*sizeof(double));
   memcpy(DTMP, DX, 6*NPAR*sizeof(double));

   if(linsolven(NPAR, 6, DTMP, MINV, A)) {
#ifdef BLAB_FIT_FAIL
      fprintf(stderr, "linsolve failure npt= %d\n", npt);
#endif
      if(TRACK_SINGLE_VEC != NULL) {
         fprintf(TRACK_SINGLE_VEC, "linsolve: %s %s %7.1f %12.5e\n",
                 vec->d[0]->id, vec->d[1]->id, (vec->d[1]->mjd-vec->d[0]->mjd)*SECDAY, sigma);
      }
      return(-1);
   }

/* NOTE: this is what is in A[]

   There are 2*NDIM components, arranged as above in S[] in x,v pairs.
   There are NPAR coefficients for each component, defined above as X[],
     nominally 0:NPAR-1 is 1,s,w,w^2,sw,s^2

   A[l+k*NPAR] is the l'th coefficient of the k'th state vector component.
 */

/* Evaluate the RMS and covariance of the fits in each coordinate */
   for(k=0; k<2*NDIM; k++) {
/* Convert DD=Sum(data*data) into an RMS wrt the fit */
      for(l=0; l<NPAR; l++) {
         DD[k] -= 2 * A[l+k*NPAR] * DX[l+k*NPAR];
         for(m=0; m<NPAR; m++)
            DD[k] += A[l+k*NPAR] * A[m+k*NPAR] * M[l+m*NPAR];
      }
      if(DD[k]>0) DD[k] = sqrt(DD[k]/npt);        // RMS, not variance
/* x,v covariance, exploiting X[0]=1 so DX[0+k*NPAR] is the sum */
      if( (k%2) == 1) {   // i.e. at k=v for x,v pair
         for(l=0; l<NPAR; l++) {
            DV[k/2] -= A[l+k*NPAR] * DX[l+(k-1)*NPAR]
               + A[l+(k-1)*NPAR] * DX[l+k*NPAR];
            for(m=0; m<NPAR; m++)
               DV[k/2] += A[l+k*NPAR] * A[m+(k-1)*NPAR] * M[l+m*NPAR];
         }
         DV[k/2] /= (npt * DD[k] * DD[k-1]);      // correlation coeff
/* But differences of small numbers might make a bogus correlation coeff... */
         if(DV[k/2] < -MAXCOV) DV[k/2] = -MAXCOV;
         if(DV[k/2] >  MAXCOV) DV[k/2] =  MAXCOV;
      }
   }

/* Save all RMS for each coordinate followed by the three covariances */
   vec->sigma = sigma;
   for(k=0; k<2*NDIM; k++) vec->rms[k] = DD[k];
   for(k=0; k<NDIM; k++)   vec->cov[k] = DV[k];

/* kd search center: average location over grid points? */
//   for(k=0; k<2*NDIM; k++)    vec->S[k] /= chinsum;
/* kd search center: infinity (s=0) position from linear fit? */
   for(k=0; k<2*NDIM; k++)    vec->S[k] = A[0+k*NPAR];

/* For each variable, save each coordinate's parameter */
   for(k=0; k<2*NDIM; k++) {
      for(l=0; l<NPAR; l++) vec->A[l+k*NPAR] = A[l+k*NPAR];
   }

   if(TRACK_SINGLE_VEC != NULL) {
      fprintf(TRACK_SINGLE_VEC, "fit for %d det with %d points: ", vec->ndet, npt);
      for(i=0; i<vec->ndet; i++) fprintf(TRACK_SINGLE_VEC, "%s%c", vec->d[i]->id,
                                         (i==vec->ndet-1)?'\n':',');
      fprintf(TRACK_SINGLE_VEC, "MJD: %10.4f sig: %10.3e S: %10.6f %10.6f %10.6f %10.3e %10.3e %10.3e\n",
              linkpar->mjd, vec->sigma, vec->S[0], vec->S[2], vec->S[4], vec->S[1], vec->S[3], vec->S[5]);
   }


/***********************************************************************/
/* If requested put out the grid with input values scattered by errors */
   if(fpscat != NULL) {
      double dr=atan(1)/45;
      double Xfit[6], chifit=0;

      fprintf(stderr, "fpscat\n");

      fprintf(fpscat, "# %s %s   unitRMS: %10.3e %10.3e %10.3e    v[dd]RMS: %10.3e %10.3e %10.3e\n",
              vec->d[0]->id, vec->d[1]->id,
              DD[0], DD[2], DD[4],
              DD[1]/dr*SECDAY, DD[3]/dr*SECDAY, DD[5]/dr*SECDAY);
      fprintf(fpscat, "# scat ir  iv    Rin      Vin   MJDeph       Reph     Veph       x[unit]      y          z        vxt[d/d]    vy        vz     Fit_x[unit]     y          z      Fit_vxt[d/d]   vy        vz       lng       lat      vlng     vlat\n");

      for(i=0; i<linkpar->nr; i++) {
         R = 1 / (smin + (i * (smax-smin)) / MAX(1,linkpar->nr-1));
         if(linkpar->nv == 0) nvgrid = i + 1;        // triangular grid
         for(j=0; j<nvgrid; j++) {
            fprintf(stderr, "fpscat %d %d\n", i, j);

            if(nvgrid == 1) {
               VR = 0.5 * (linkpar->vmin + linkpar->vmax);
            } else {
               VR = linkpar->vmin + (j * (linkpar->vmax-linkpar->vmin)) / (nvgrid-1);
            }
            R -= VR*1e3 * Dt / AU;
            if(R < RMIN) R = RMIN;   // Floor for R

            lnr = log(R);
            wr = VR*1e3/(R*AU)*SECDAY;

            chifit += puma_chi(puma, lnr, wr, RMS);
            puma_bary(puma, ephdet, linkpar->integrate, &mjd, X, V, dX, dV, &sigma);
            rho = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
            rhodot = (X[0]*V[0]+X[1]*V[1]+X[2]*V[2]) / rho;

            if(ABS(rhodot/rho*Dt) > linkpar->wdtmin) continue;
            if( (v2=(V[0]*V[0]+V[1]*V[1]+V[2]*V[2])) > linkpar->vtmax*linkpar->vtmax*1e6) continue;

/* Evaluate the fit for this r,vr: S[] = three components of x/r,vx/r */
            s = 1 / (rho/AU);             // [AU^-1]     1/r
            w = (rhodot/1e3) / (rho/AU);  // [km/s/AU]  vr/r
            Xfit[0] = 1;
            Xfit[1] = s;
            Xfit[2] = w;
            Xfit[3] = w*w;
            Xfit[4] = s*w;
            Xfit[5] = s*s;
            for(k=0; k<2*NDIM; k++) {
               for(l=0, S[k]=0.0; l<NPAR; l++) S[k] += A[l+k*NPAR] * Xfit[l];
            }

/* Make S[] v tangential only */
            S[1] -= S[0] * rhodot/rho;
            S[3] -= S[2] * rhodot/rho;
            S[5] -= S[4] * rhodot/rho;
//            printf("%8.4f %8.2f %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", s, w, S[0], S[1], S[2], S[3], S[4], S[5]);

            double dlnr=0.0, dwr=0.0, slnr, swr;
            int iscat;
            for(iscat=0; iscat<nscat; iscat++) {
/* Angular versions of unit and velocity vectors */
               double TH[NDIM], PH[NDIM], rxy;
               double vlng, vlat, lng, lat;
               lng = atan2(X[1],X[0]);
               if(lng < 0) lng += 8*atan(1);
               lat = asin(X[2]/rho);
               rxy = sqrt(X[0]*X[0]+X[1]*X[1]);
               TH[0] = X[0]*X[2] / rxy / rho;
               TH[1] = X[1]*X[2] / rxy / rho;
               TH[2] = -rxy / rho;                           // wrt pole
               PH[0] = -X[1] / rxy;
               PH[1] =  X[0] / rxy;
               PH[2] = 0;
               vlng =  (PH[0]*V[0]+PH[1]*V[1]+PH[2]*V[2]) / rho;
               vlat = -(TH[0]*V[0]+TH[1]*V[1]+TH[2]*V[2]) / rho; // lat

/* Make V tangential only */
               V[0] -= X[0]/rho * rhodot;
               V[1] -= X[1]/rho * rhodot;
               V[2] -= X[2]/rho * rhodot;
               
               fprintf(fpscat, "%4d %4d %3d %8.6f %7.2f  %7.1f    %8.5f %8.4f   %10.7f %10.7f %10.7f  %9.6f %9.6f %9.6f  %10.7f %10.7f %10.7f  %9.6f %9.6f %9.6f  %9.5f %9.5f %8.5f %8.5f\n", 
                       iscat, i, j+1, R, VR, mjd,      // initial R,VR
                       rho/AU, rhodot/1e3,             // [AU, km/s] at mjd
                       X[0]/rho, X[1]/rho, X[2]/rho,   // [] unit vector
                       V[0]/rho/dr*SECDAY,             // [deg/day] tangvel
                       V[1]/rho/dr*SECDAY, V[2]/rho/dr*SECDAY,
                       S[0], S[2], S[4],               // [] fit unit vector
                       S[1]/dr*SECDAY,                 // [deg/day] fit tangvel
                       S[3]/dr*SECDAY, S[5]/dr*SECDAY,
                       lng/dr, lat/dr,                 // [deg] lng,lat
                       vlng/dr*SECDAY, vlat/dr*SECDAY); // [deg/day] tangvel
/* reload X,V with new values by perturbing input data */
               slnr = lnr;
               swr = wr;
               chin = puma_perturb(puma, &slnr, &swr, dlnr, dwr, RMS);
               puma_bary(puma, ephdet, linkpar->integrate, &mjd, X, V, dX, dV, &sigma);
               rho = sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
               rhodot = (X[0]*V[0]+X[1]*V[1]+X[2]*V[2]) / rho;
            }
         }
      }
   }
/* end of scatter output */
/***********************************************************************/




#ifdef RESULTS_TO_STDOUT_FOR_DEBUGGING
/* Start output: vec->d[0] vec->d[1] dt Dt sigma */
   printf("%s %s %7.1f %8.1f %12.5e",
          vec->d[0]->id, vec->d[1]->id, 
          (vec->d[1]->mjd-vec->d[0]->mjd)*SECDAY,
          (linkpar->mjd-0.5*(vec->d[0]->mjd+vec->d[1]->mjd))*SECDAY,
          sigma);

/* Write all RMS for each coordinate followed by the three covariances */
   for(k=0; k<2*NDIM; k++) printf(" %10.3e", DD[k]);
   for(k=0; k<NDIM; k++) printf(" %6.3f", DV[k]);

/* For each variable, write each coordinates parameter */
   for(l=0; l<NPAR; l++) {
      for(k=0; k<2*NDIM; k++) {
         if(l==0) printf(" %12.5e", A[0+k*NPAR]);
         else     printf(" %10.3e", A[l+k*NPAR]);
      }
   }
   printf("\n");
#endif


   return(0);
}

/* Solve nxn matrix equation Y = M X, return X for ny versions (mangle M,Y!) */
int linsolven(int n, int ny, double *Y, double *M, double *X)
{
   int i, j, k, iy, *rowstatus, *row;
   double rat;

   rowstatus = (int *)calloc(n, sizeof(int));
   row = (int *)calloc(n, sizeof(int));

/* Solve matrix by Gaussian elimination */
   for(j=0; j<n; j++) {		/* j = column to be zero'ed out */

/* Choose a good equation to work on (pivot row) */
      for(i=k=0, rat=0.0; i<n; i++) {
	 if(rowstatus[i]) continue;
	 if(ABS(M[j+i*n]) > rat) {
	    k = i;
	    rat = ABS(M[j+k*n]);
	 }
      }
      rowstatus[k] = 1;
      row[j] = k;
      if(rat == 0.0) {  // singular
         free(row);
         free(rowstatus);
	 return(-1);
      }
/* Subtract away matrix below diagonal */
      for(i=0; i<n; i++) {	/* i = subsequent equations */
	 if(rowstatus[i]) continue;
	 rat = M[j+i*n] / M[j+row[j]*n];
	 for(iy=0; iy<ny; iy++) Y[i+iy*n] -= rat * Y[row[j]+iy*n];
 	 for(k=j+1; k<n; k++) M[k+i*n] -= rat * M[k+row[j]*n];
      }
   }

/* Back substitute into upper diagonal matrix to solve for parameters */
   for(j=n-1; j>=0; j--) {
      for(iy=0; iy<ny; iy++) {
	 X[j+iy*n] = Y[row[j]+iy*n];
	 for(k=j+1; k<n; k++) X[j+iy*n] -= M[k+row[j]*n] * X[k+iy*n];
	 X[j+iy*n] /= M[j+row[j]*n];
      }
   }
   free(row);
   free(rowstatus);
   return(0);
}

/* Angle betweeen two points on the sphere, all [deg] */
double sphdist(double r1, double d1, double r2, double d2)
{
   double x1, x2, y1, y2, z1, z2, dot, dr=atan(1)/45;
   x1 = cos(d1*dr) * cos(r1*dr);
   y1 = cos(d1*dr) * sin(r1*dr);
   z1 = sin(d1*dr);
   x2 = cos(d2*dr) * cos(r2*dr);
   y2 = cos(d2*dr) * sin(r2*dr);
   z2 = sin(d2*dr);
   dot = x1*x2 + y1*y2 + z1*z2;
   if(dot > 1) return(0.0);
   if(dot < -1) return(180.0);
   return(acos(dot)/dr);
}

#define VECALLOC 1000

/* Detection pairs are guaranteed to be successive in time */
/* Identify pairs of detections within ang vel w[deg/day] of each other */
int link_det(LINKPAR *linkpar, int ndet, DET *det, int epoch, VEC **vec)
{
   int j, nvec, nalloc;
   void *kd, *set;
   double P1[4], P2[4], d2, du, dth, dr=atan(1)/45;
   DET *dp2;

/* Initialize vector count */
   *vec = NULL;
   nvec = nalloc = 0;

/* Create a 4D structure for the unit vectors and observation time */
   kd = kd_create(4);
   
/* Pair the detections within a given angular velocity */
   du = linkpar->omega * dr * linkpar->dtmax;    // [rad] max dist to search

/* Insert all the unit vectors into a kd-tree */
   for(j=0; j<ndet; j++) {

      if(TRACK_SINGLE_PAIR != NULL) {
         fprintf(TRACK_SINGLE_PAIR, "det insert %s\n", det[j].id);
      }

/* Just from epoch 1 or epoch 2? */
      if(epoch == 1 && linkpar->split && det[j].mjd > linkpar->mjd) continue;
      if(epoch == 2 && linkpar->split && det[j].mjd <= linkpar->mjd) continue;

      P2[0] = cos(det[j].dec*dr) * cos(det[j].ra*dr);
      P2[1] = cos(det[j].dec*dr) * sin(det[j].ra*dr);
      P2[2] = sin(det[j].dec*dr);
      P2[3] = (det[j].mjd - linkpar->mjd) / linkpar->dtmax * du;
      kd_insert(kd, P2, (void *)&det[j]);
   }

/* Loop over all elements looking for proximity */
   for(j=0; j<ndet; j++) {

      if(TRACK_SINGLE_PAIR != NULL) {
         fprintf(TRACK_SINGLE_PAIR, "det pair %s\n", det[j].id);
      }

/* Just from epoch 1 or epoch 2? */
      if(epoch == 1 && linkpar->split && det[j].mjd > linkpar->mjd) continue;
      if(epoch == 2 && linkpar->split && det[j].mjd <= linkpar->mjd) continue;

      P1[0] = cos(det[j].dec*dr) * cos(det[j].ra*dr);
      P1[1] = cos(det[j].dec*dr) * sin(det[j].ra*dr);
      P1[2] = sin(det[j].dec*dr);
      P1[3] = (det[j].mjd - linkpar->mjd) / linkpar->dtmax * du;

/* Find detections within range of this detection */
      set = kd_nearest_range(kd, P1, du);

/* Examine each one, testing for distance less than w*|mjd1-mjd2| */
      while( !kd_res_end(set) ) {

/* get the data and position of the current result item */
         dp2 = (DET *)kd_res_item(set, P2);

/* Skip dups or links prior to det[j] */
         if(dp2 <= det+j) {
            kd_res_next(set);   // next entry
            if(TRACK_SINGLE_PAIR != NULL) {
               fprintf(TRACK_SINGLE_PAIR, "det is prior %s\n", det[j].id);
            }

            continue;
         }

/* Actually allowable dtheta for this actual difference in MJD */
         dth = ABS(dp2->mjd-det[j].mjd) * linkpar->omega * dr;

/* Close enough?  Save it. */
         d2 = (P1[0]-P2[0])*(P1[0]-P2[0]) +
              (P1[1]-P2[1])*(P1[1]-P2[1]) +
              (P1[2]-P2[2])*(P1[2]-P2[2]);
         if(TRACK_SINGLE_PAIR != NULL) {
            fprintf(TRACK_SINGLE_PAIR, "det later %s %10.3e %10.3e %s\n", det[j].id, dth, sqrt(d2), dp2->id);
         }
         if(d2 < dth*dth) {
            if(nvec >= nalloc) {
               nalloc += VECALLOC;
               *vec = (VEC *)realloc(*vec, nalloc*sizeof(VEC));
            }
            (*vec)[nvec].d = (DET **)calloc(2, sizeof(DET *));
            (*vec)[nvec].ndet = 2;
            (*vec)[nvec].d[0] = det + j;
            (*vec)[nvec].d[1] = dp2;
            (*vec)[nvec].mjd = 0.5*(det[j].mjd + dp2->mjd);
            (*vec)[nvec].dt = ABS(dp2->mjd - det[j].mjd) * SECDAY;
            nvec += 1;
         }

         kd_res_next(set);   // next entry
      }
      kd_res_free(set);
   }

   kd_free(kd);

   return(nvec);
}

/* Look up a string among a sorted lists of strings, return index */
int string_lookup(char *s, int n, char **list)
{
   int i, i0=0, i1=n-1, g;

/* Not there or matches end? */
   g = strcmp(s, list[0]);
   if(g < 0)  return(n);
   if(g == 0) return(0);
   g = strcmp(s, list[n-1]);
   if(g > 0)  return(n);
   if(g == 0) return(n-1);

/* Binary chop, i0 - i1, not on end */
   while(i1 > i0+1) {
      i = (i0+i1) / 2;
      g = strcmp(s, list[i]);
      if(g == 0) return(i);
      if(g < 0) i1 = i;
      else      i0 = i;
   }
   return(n);
}



/* Read an input list of detections and create vectors */
/* Potentially augments a vector list from link_det() */
int read_vec(char *infile, LINKPAR *linkpar, int ndet, DET *det,
             int nvec, VEC **vec)
{
   int i, j, n, nv0, nalloc, nd;
   char *csv, *id, *c;
   double mjd, m1=0, m2=0;
   FILE *fp;

/* Open input */
   if(strcmp(infile, "-") == 0) {
      fp = stdin;
   } else if( (fp = fopen(infile, "r")) == NULL) {
      fprintf(stderr, "Cannot open %s for reading\n", infile);
      exit(1);
   }

   nv0 = nalloc = nvec;
   csv = (char *)malloc(102400);
   id = (char *)malloc(102400);

/* Sort all the detections by ID for lookup of vector references */
   char **srtid = calloc(ndet, sizeof(char *));
   int *srtdet = calloc(ndet, sizeof(int));
   for(i=0; i<ndet; i++) {
      srtid[i] = det[i].id;
      srtdet[i] = i;
   }
   tsort_str(ndet, srtid, srtdet);


/* Read all the data */
   for(n=0;   ; n++) {
      if(fgets(csv, 102400, fp) == NULL) break;
      if(csv[0] == '#') continue;                       // Skip comments
      if(strlen(csv) == 0) continue;                    // Skip empty lines

/* How many commas and detections in this line? */
      for(c=csv, nd=1; c!=NULL && c<csv+strlen(csv)-2; nd++) {
         if( (c=index(c+1, ',')) == NULL) break;
         c++;
      }
      if(nd < 2) {
         fprintf(stderr, "Vector %s need at least 2 detections\n", csv);
         exit(1);
      }
//      fprintf(stderr, "csv= %s nd= %d\n", csv, nd);

/* Make space */
      if(nvec >= nalloc) {
         nalloc = nvec + VECALLOC;
         *vec = (VEC *)realloc(*vec, nalloc*sizeof(VEC));
      }
/* Save the detection pointers */
      (*vec)[nvec].ndet = nd;
      (*vec)[nvec].d = (DET **)calloc(nd, sizeof(DET *));
      mjd = 0;
/* Find each ID from the list amidst all the detections */
      for(c=csv, i=0; i<nd; i++) {
         sscanf(c, "%[^,\n]", id);
//         fprintf(stderr, "id= %s\n", id);

#ifdef SLOW_BUT_COMPACT
/* Look up this ID among all the detections */
         for(j=0; j<ndet && strcmp(id, det[j].id)!=0; j++);
#else
/* Look up this ID among all the detections, dereference if legal */
         if( (j = string_lookup(id, ndet, srtid)) < ndet) j = srtdet[j];
#endif
         if(j == ndet) {
            fprintf(stderr, "No detection for vector id %s from %s", id, csv);
            exit(1);
         }
         (*vec)[nvec].d[i] = det + j;
         mjd += det[j].mjd;
         if(i == 0) {
            m1 = m2 = det[j].mjd;
         } else {
            m1 = MIN(m1, det[j].mjd);
            m2 = MAX(m2, det[j].mjd);
         }

         c = index(c, ',') + 1;
      }
      (*vec)[nvec].mjd = mjd / nd;
      (*vec)[nvec].dt = (m2 - m1) * SECDAY;
//      fprintf(stderr, "vector %d mjd %.1f dt %.1f\n", nvec, (*vec)[nvec].mjd, (*vec)[nvec].dt);

      nvec++;
   }

   if(strcmp(infile, "-") != 0) fclose(fp);

   free(srtid);
   free(srtdet);

   free(csv);
   free(id);

   return(nvec-nv0);
}





#define DETALLOC 1000

/* read_trd() reads a TRD9 file, sorts these detections by MJD (primary)
 * and id (secondary) and saves the results in a DET array
 */

/* Read an input file of trd data */
int read_trd(char *infile, DET **det, LINKPAR *linkpar)
{
   int i, i1, i2, n, ndet, nitem, nalloc;
   char line[1024], id[1024];
   DET d1, *dread;
   FILE *fp;

/* Open input */
   if(strcmp(infile, "-") == 0) {
      fp = stdin;
   } else if( (fp = fopen(infile, "r")) == NULL) {
      fprintf(stderr, "Cannot open %s for reading\n", infile);
      exit(1);
   }

   ndet = nalloc = 0;
   dread = NULL;

/* Read all the data */
   for(n=0;   ; n++) {
      if(fgets(line, 1024, fp) == NULL) break;
      if(line[0] == '#') continue;	/* Skip comments */
/* Read a trd */
      nitem = sscanf(line, "%lf %lf %lf %f %f %f %f %f %s",
                     &d1.mjd, &d1.ra, &d1.dec, &d1.xerr, &d1.terr,
                     &d1.lng, &d1.lat, &d1.elev, id);

      if(nitem <= 0) continue;
      if(nitem < 9) {
	 fprintf(stderr, "Error reading from `%s' at line %d\n", infile, n);
	 fprintf(stderr, "Line `%s'\n", line);
	 exit(1);
      }

/* Add d1 to the detection array */
      if(ndet >= nalloc) {
         nalloc = ndet + DETALLOC;
         dread = (DET *)realloc(dread, nalloc*sizeof(DET));
      }
      dread[ndet] = d1;
      dread[ndet].id = malloc(strlen(id)+1);
      strcpy(dread[ndet].id, id);
      dread[ndet].dgrp = -1;
      if(TRACK_SINGLE_DET != NULL) {
         fprintf(TRACK_SINGLE_DET, "det %s %d\n", dread[ndet].id, ndet);
      }
      ndet += 1;
   }
   if(strcmp(infile, "-") != 0) fclose(fp);

/* Sort all the detections by MJD (primary) and ID (secondary) */
   double *mjdsrt = (double *)calloc(ndet, sizeof(double));
   int *idx = (int *)calloc(ndet, sizeof(int));
   for(i=0; i<ndet; i++) {
      mjdsrt[i] = dread[i].mjd;
      idx[i] = i;
   }
   tsort_d(ndet, mjdsrt, idx);

   *det = (DET *)calloc(ndet, sizeof(DET));
   nalloc = DETALLOC;
   char **csrt = (char **)calloc(nalloc, sizeof(char *));

/* Squeeze out all the dups */
   for(i1=0; i1<ndet; i1=i2) {
      for(i2=i1+1; i2<ndet && dread[idx[i1]].mjd==dread[idx[i2]].mjd; i2++);

      if(i2 == i1+1) {
         (*det)[i1] = dread[idx[i1]];

      } else {
         if(i2-i1 >= nalloc) {
            nalloc = i2-i1 + DETALLOC;
            csrt = (char **)realloc(csrt, nalloc*sizeof(char *));
         }
         for(i=0; i<i2-i1; i++) csrt[i] = dread[idx[i+i1]].id;
         tsort_str(i2-i1, csrt, idx+i1);
         for(i=0; i<i2-i1; i++) (*det)[i+i1] = dread[idx[i+i1]];
      }
   }
/* Assign a unique index to each detection */
   for(i=0; i<ndet; i++) (*det)[i].idx = i;
   linkpar->ndet = ndet;

   if(linkpar->mjd == 0) {
/* Find all groups that link by dtmax and the biggest gap between them */
      double mjd1, gap0, gap1;        // group and gap delimiters
      mjd1 = gap0 = gap1 = mjdsrt[0];
      for(i=1; i<ndet; i++) {
         if(mjdsrt[i] <= mjd1+linkpar->dtmax && i<ndet-1) {
            mjd1 = mjdsrt[i];
         } else {
//            printf("%9.3f %9.3f   %9.3f %9.3f", mjd0, mjd1, gap0, gap1);
            if(mjdsrt[i]-mjd1 > gap1-gap0) {
               gap0 = mjd1;
               gap1 = mjdsrt[i];
            }
            mjd1 = mjdsrt[i];
//            printf(" -> %9.3f %9.3f\n", gap0, gap1);
         }
      }
/* Set the reference MJD to half way through the biggest gap */
      linkpar->mjd = 0.5 * (gap0 + gap1);
      linkpar->mjd = linkpar->mfrac * NINT(linkpar->mjd/linkpar->mfrac);
   }

   free(dread);
   free(mjdsrt);
   free(idx);
   free(csrt);
   return(ndet);
}
