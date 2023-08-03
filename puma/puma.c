/* Fit a solar system orbit to data: Position Using Motion with Acceleration */
/* v1.21 - 200101 separated into library and application */
/* v1.2  - 191212 modified to do a much better job on longer arcs */
/* v1.0  - 180605 John Tonry, simplified from ssorb.c */


/* Some examples...

  cd /atlas/catdev/PUMA
  rock=2019WR3
  puma2 -in $rock.trd -r 0.04906021 -v -2.961 -resid -


  for rock in $(basename -a -s .trd *.trd) ; do 
    echo ; echo $rock "chi/N"
    puma2 -in $rock.trd -fit -resid -
  done | grep chi/N
*/


/* Get examples from MPC:

  https://minorplanetcenter.net/db_search

e.g. the following have multiple reports from ATLAS

  2019 WR3
  2019 WQ3
  2019 XC1
  2019 WR2
  2006 WH1

Unpack MPC, ugh

  rock=2019XC1

  awk '{date=substr($0,16,17);
	Y=substr(date,1,4); M=substr(date,6,2); D=substr(date,9,2);
	h=int(24*substr(date,11)); m=int(60*(24*substr(date,11)-h)); 
	s=60*(60*(24*substr(date,11)-h)-m); 
	printf "%04d-%02d-%02dT%02d:%02d:%05.2f\n", Y,M,D,h,m,s}' 2019XC1.mpc | date -u -f - "+%s.%N" | awk '{printf "%.7f\n",$1/86400.0+40587.0}' > /tmp/$rock.mjd

  awk '{mpnum=substr($0,1,5); name=substr($0,6,7); date=substr($0,16,17);
        ra=substr($0,33,12); dsgn=substr($0,45,1); dec=substr($0,46,12); 
        m=substr($0,66,5); filt=substr($0,71,1); obs=substr($0,78,3);
	r=15*(substr(ra,1,2)+substr(ra,4,2)/60+substr(ra,7)/3600);
	d=(substr(dec,1,2)+substr(dec,4,2)/60+substr(dec,7)/3600);
	if(dsgn == "-") d=-d;
	if(obs == "T05") {lng=-156.2570; lat=20.7076; alt=3041}
	if(obs == "T08") {lng=-155.5763; lat=19.5363; alt=3412}
	err=0.1;
	printf "%9.5f %9.5f %5.2f %5.2f %9.4f %9.4f %4d %s %5.2f %s %s\n", 
             r, d, err, err, lng, lat, alt, name, m, filt, obs}' $rock.mpc |
	paste -d ' ' /tmp/$rock.mjd - > $rock.trd

  rm /tmp/$rock.mjd
*/

/* TODO:
 o why does phi(t) not have a good zero?
 * iteration to correct rho is just not right...
 * 2019WR3 great circle residuals?!  Needs a cubic, implement poly_gctfit()
 * modify grid to respect vr/r: [deg/day]
 * return/print tangential velocity?
 * closest approach time, pos, etc.
 * great circle residuals in both resid tables
 * library
 * -ephmjd MJD1:MJD2:dtime  add to man page
 * r-vr covariance!  what's going on???
 * separate det structures for data and ephemerides
 * separate resid tables for data and ephemerides: resid table different
 * acc towards barycenter (or include moon)
 * move to barycenter
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "puma.h"

// #define UNWEIGHTED_XROCK	/* Define to use no weights in trajectory fit */

#define MAXEPHMJD 5000	/* Max number of ephemeris dates */
#define SECDAY (86400)		/* Seconds in a day */
#define AU (1.49598E11) 	/* [m] */
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/* Global variables */
extern int PUMA_VERB;

int JSON;       	/* JSON format for output */


/* Read an input file, return number of detections and data */
int rock_input(char *infile,
               double lng_default, double lat_default, double elev_default,
               double lng_eph, double lat_eph, double elev_eph,
               double poserr_default, char *epharg, void *p);

/* Parse a CSV string with real variables, return array and number */
int parse_eph(char *arg, int maxvar, double *var, char **id,
              double *lng, double *lat, double *elev);

/* Uncertainty in Geo E error, given lnr and wr=vr/r */
double geo_E_error(double lnr, double dlnr, double wr, double dwr,
                   double rwcov, double E);

/* Quicksort doubles */
void tsort_d(int n, double *x, int *idx);


void usage(int error) 
{
/* Emit usage message */
   fprintf(stderr, \
"Usage: puma -in INFILE [-fit] [-grid] [args] [-verbose LEVEL]\n\
\n\
  -in INFILE : input file, can be - for STDIN\n\
  -resid fname : write residuals for each detection to fname\n\
  -fit : non-linear least squares fit to detections\n\
  -grid : grid search for best fit\n\
  -gridout fname : write grid search results to fname (or 'json')\n\
\n  -ephmjd M[,M2,M3,...] : integrate trajectory to MJD M\n\
  -eph fname : scatter plot file name\n\
  -json : write json format\n\
  \n-verbose LEVEL : set debug verbosity to LEVEL (0,1,2)\n\
\n\
`man puma' for full documentation.\n\
");
   exit(error);
}


extern int NXROCK;		// count of calls to xrock()
extern int DEBUG;


int main(int argc, char **argv)
{
   int i, j, k, npt, nrad=0, NephMC, dofit, eph, wgrid, geobound;
   char *infile, *ephout, *resout, *scatout, *gridout, *epharg;
   char line[1024];
   double R, VR, dt, chi, chi0, chimax;
   double *chisrt=NULL, *chirad=NULL, *chivel=NULL;
   double obslng, obslat, obselev, obserr, ephlng, ephlat, ephelev;
   int *chidx=NULL;

   FILE *fpeph, *fpres, *fpgrid, *fpscat;
   int nrgrid, nvgrid;
   double rmin, rmax, vmin, vmax, wmin, wmax, RFIT, VFIT;
   
   void *puma;
   double t0, lnr, wr, dlnr, dwr, rwcov, raobs, decobs, mjd, rms[4];
   double raperi, decperi, mjdperi, rperi, wperi, ragct, decgct;
   double sigma, rafit, decfit, delta, vrad, wtan, vpa, xerr, terr, slnr, swr;
   double Esun, Egeo, dEgeo, E, ecc, incl, a, period;
//   double vsun, vrad0, vtan0;

   int nfit, ndet, whichorigin, ndate;
   double X[3], V[3], A[9], gcpole[3], gcorigin[3], gcphi[PUMA_MAXPHI];

   struct timeval tv0, tv1;

/* Start time */
   gettimeofday(&tv0, NULL);


   int iter, niter;

///////////////////////////
///      ARGUMENTS      ///
///////////////////////////
/* Defaults */
   infile = "-";	/* Input data filename */
   resout = NULL;	/* List of residuals */
   ephout = NULL;	/* Ephemerides output file */
   gridout = NULL;	/* Grid output file */
   scatout = NULL;	/* Statistics scatter output file */
   fpres = NULL;	/* Residual file */
   fpeph = NULL;	/* Ephemerides file */
   fpgrid = NULL;	/* Grid file */
   fpscat = NULL;	/* Statistics scatter file */
   R = -1;		/* [m] obs distance for linearized fit */
   VR = 0.0;		/* [m/s] obs radial velocity linearized fit */
   dofit = 0;		/* Do a fit instead of evaluation at specific R,V */
   NephMC = 0;		/* Compute ephemerides uncertainties? = N random pts */
   obslng = obslat = obselev = -1000.0;	/* Observer location */
   obserr = 0.1;	/* [arcsec] nominal observational error */
   ephlng = ephlat = ephelev = -1000.0;	/* Ephemerides location */
   epharg = NULL;	/* ephemerides input argument */
   chimax = 20;		/* max chi^2/N for scatter file point */

   nrgrid = 0;		/* Number of grid search in r */
   nvgrid = 0;		/* Number of grid search in vr */
   wgrid = 0;		/* Grid in w or v? */
   rmin = 1e-4;		/* [AU] Minimum r */
   rmax = 10;		/* [AU] Maximum r */
   vmin = -20;		/* [km/s] Minimum vr */
   vmax = 20;		/* [km/s] Maximum vr */

   niter = 0;		/* iterations */


/* Parse the arguments */
   if (argc < 2) {
      usage(1);         /* usage then exit error */
   }
   for(i=1; i<argc; i++) {
      if (strcmp(argv[i], "-verbose") == 0) {	   /* debug? */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%d", &PUMA_VERB);

      } else if (strcmp(argv[i], "-in") == 0) {	   /* input file */
         if (i == argc - 1) usage(1);
	 infile = argv[++i];

      } else if (strcmp(argv[i], "-lng") == 0) {   /* [deg] observatory long */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%lf", &obslng);

      } else if (strcmp(argv[i], "-lat") == 0) {   /* [deg] observatory lat */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%lf", &obslat);

      } else if (strcmp(argv[i], "-elev") == 0) {   /* [m] observatory elev */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%lf", &obselev);

      } else if (strcmp(argv[i], "-ephlng") == 0) {   /* [deg] ephem long */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%lf", &ephlng);

      } else if (strcmp(argv[i], "-ephlat") == 0) {   /* [deg] ephem lat */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%lf", &ephlat);

      } else if (strcmp(argv[i], "-ephelev") == 0) {   /* [m] ephem elev */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%lf", &ephelev);

      } else if (strcmp(argv[i], "-obserr") == 0) { /* [arcsec] posn error */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%lf", &obserr);

      } else if (strcmp(argv[i], "-json") == 0) { /* JSON output? */
	 JSON = 1;

      } else if (strcmp(argv[i], "-fit") == 0) {  /* Do a fit? */
	 dofit = 1;

      } else if (strcmp(argv[i], "-r") == 0) {  /* [AU] linear initial dist */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%lf", &rmin);
	 rmax = rmin;
	 nrgrid = 1;

      } else if (strcmp(argv[i], "-v") == 0) {  /* [km/s] linear initial vr */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%lf", &vmin);
	 vmax = vmin;
	 nvgrid = 1;

      } else if (strcmp(argv[i], "-grid") == 0) {  /* grid spec */
	 nrgrid = nvgrid = 20;		/* Grid count for nominal limits above */
	 if(i<argc-1 && sscanf(argv[i+1], "%d,%lf,%lf,%d,%lf,%lf",
		       &nrgrid, &rmin, &rmax, &nvgrid, &vmin, &vmax) == 6) {
	    i++;
	 }

      } else if (strcmp(argv[i], "-wgrid") == 0) {  /* grid in vr/r? */
	 wgrid = 1;

      } else if (strcmp(argv[i], "-gridout") == 0) { /* Grid output file */
         if (i == argc - 1) usage(1);
	 gridout = argv[++i];

      } else if (strcmp(argv[i], "-resid") == 0) { /* Residual output file */
         if (i == argc - 1) usage(1);
	 resout = argv[++i];

      } else if (strcmp(argv[i], "-ephmjd") == 0) {/* ephem at particular MJDs */
         if (i == argc - 1) usage(1);
	 epharg = argv[++i];

      } else if (strcmp(argv[i], "-eph") == 0) {   /* Ephemerides output file */
         if (i == argc - 1) usage(1);
	 ephout = argv[++i];

      } else if (strcmp(argv[i], "-scat") == 0) {   /* Stat scatter output file */
         if (i == argc - 1) usage(1);
	 scatout = argv[++i];

      } else if (strcmp(argv[i], "-nscat") == 0) {/* number of random eph err */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%d", &NephMC);

      } else if (strcmp(argv[i], "-chimax") == 0) {/* max chi for scatter */
         if (i == argc - 1) usage(1);
	 sscanf(argv[++i], "%lf", &chimax);

      } else if (strcmp(argv[i], "-ndate") == 0) {/* acc(t) fit order */
//	 sscanf(argv[++i], "%d", &NDATE);

      } else if (strcmp(argv[i], "-niter") == 0) {   /* timing iterations */
	 sscanf(argv[++i], "%d", &niter);

      } else {
	 fprintf(stderr, "Unknown arg `%s'\n", argv[i]);
	 exit(1);
      }
   }

/* If a scatter file, use a default for NephMC */
   if(scatout != NULL && NephMC == 0) {
      NephMC = 100;
   }

/* If NephMC, default for scatter file is stdout */
   if(scatout == NULL && NephMC > 0) {
      fpscat = stdout;		// Default output is stdout
   }

/* If NephMC we must have a fit */
   if(NephMC) dofit = 1;

/* If NephMC we must have a target date */
   if(NephMC && epharg == NULL) {
      fprintf(stderr, "puma error: ephemeris error enabled without a target -ephmjd date\n");
      exit(1);
   }

/* If no fit and no grid requested, use a 3x3 grid as a default */
   if(nrgrid == 0 && !dofit) nrgrid = nvgrid = 3;


///////////////////////////
///      ARGUMENTS      ///
///////////////////////////
   


/////////////////////
///      I/O      ///
/////////////////////
/* Open grid file if requested */
   if(gridout != NULL &&
      (strcmp(gridout, "-") == 0 || 
       strcmp(gridout, "json") == 0)) {
      fpgrid = stdout;
   } else {
      if(gridout != NULL && (fpgrid = fopen(gridout, "w")) == NULL) {
	 fprintf(stderr, "Cannot open %s for writing\n", gridout);
	 exit(1);
      }
   }

/* Open residual file if requested */
   if(resout != NULL && strcmp(resout, "-") == 0) {
      fpres = stdout;
   } else {
      if(resout != NULL && (fpres = fopen(resout, "w")) == NULL) {
	 fprintf(stderr, "Cannot open %s for writing\n", resout);
	 exit(1);
      }
   }

/* Open output ephemeris best fit file if requested */
   if(ephout != NULL && (
	 strcmp(ephout, "-") == 0 || strcmp(ephout, "json") == 0)) {
      fpeph = stdout;
   } else {
      if(ephout != NULL && (fpeph = fopen(ephout, "w")) == NULL) {
	 fprintf(stderr, "Cannot open %s for writing\n", ephout);
	 exit(1);
      }
   }

/* Open output ephemeris scatter file if requested */
   if(scatout != NULL && (
	 strcmp(scatout, "-") == 0 || strcmp(scatout, "json") == 0)) {
      fpscat = stdout;
   } else {
      if(scatout != NULL && (fpscat = fopen(scatout, "w")) == NULL) {
	 fprintf(stderr, "Cannot open %s for writing\n", scatout);
	 exit(1);
      }
   }

/////////////////////
///      I/O      ///
/////////////////////



////////////////////////////
///      INITIALIZE      ///
////////////////////////////
   if(puma_init(&puma)) {
      fprintf(stderr, "Cannot initialize puma\n");
      exit(1);
   }



////////////////////////////
///      INPUT DATA      ///
////////////////////////////
/* Read the detection data */
   npt = rock_input(infile,
                    obslng, obslat, obselev, ephlng, ephlat, ephelev,
                    obserr, epharg, puma);

/* Get the number of fit points */
   puma_fitpar(puma, &nfit, &ndet, &whichorigin, &ndate,
	       &t0, X, V, A, gcpole, gcorigin, gcphi);

   if(nfit < 2 && (nrgrid != 1 || nvgrid != 1) ) {
      fprintf(stderr, "Cannot evaluate an orbit with only %d points\n", nfit);
      exit(1);
   }
   if(nfit < 3 && dofit) {
      fprintf(stderr, "Cannot fit an orbit with only %d points\n", nfit);
      exit(1);
   }
   if(PUMA_VERB > 0) printf("Read %d observations and %d ephemeris dates\n",
		       nfit, npt-nfit);

/* Ephemeris dates imply write to stdout if not otherwise specified */
   if(npt > nfit && fpeph == NULL) fpeph = stdout;


////////////////////////////
///      INPUT DATA      ///
////////////////////////////

   if(JSON) {
      printf("{ \n");
   }

////////////////////////////////////////////////////////////////
/* Iterations to test timing */
   for(iter=0; iter<MAX(1,niter); iter++) {
////////////////////////////////////////////////////////////////

/* Defaults for fit */
   RFIT = sqrt(rmin*rmax);
   VFIT = 0.5*(vmin+vmax) + 1e-2;
//   printf("%10.3e %10.3e %10.3e  %6.2f %6.2f %6.2f\n",
//	  rmin, rmax, RFIT, vmin, vmax, VFIT);



/////////////////////////////
///      GRID SEARCH      ///
/////////////////////////////
   if(nrgrid >= 1 || nvgrid >= 1) {

      if(JSON) {
	 printf("\"nrgrid\": %d, \"nvgrid\": %d, \"grid\": [\n", nrgrid, nvgrid);
      } else if(fpgrid != NULL) {
	 fprintf(fpgrid, "# ir iv    R[AU]     logR[AU] VR[km/s] VR/R[1/d]   chi/N\n");
      }		       

      chi0 = 1e38;
      chisrt = (double *)calloc(nrgrid*nvgrid, sizeof(double));
      chidx = (int *)calloc(nrgrid*nvgrid, sizeof(int));
      chirad = (double *)calloc(nrgrid*nvgrid, sizeof(double));
      chivel = (double *)calloc(nrgrid*nvgrid, sizeof(double));

/* Step over the entire grid, remembering best fit */
      wmin = (vmin<0) ? vmin/rmin : vmin/rmax;
      wmax = (vmax>0) ? vmax/rmin : vmax/rmax;
      for(i=0; i<nrgrid; i++) {
	 for(j=0; j<nvgrid; j++) {
	    R = rmin * exp((i*log(rmax/rmin))/MAX(1,nrgrid-1));	// [AU]
	    if(wgrid) {
	       VR = R * (wmin + (j*(wmax-wmin))/MAX(1,nvgrid-1));// [km/s]
	    } else {
	       VR = vmin + (j*(vmax-vmin))/MAX(1,nvgrid-1);	// [km/s]
	    }

	    lnr = log(R);
	    wr = VR*1e3/(R*AU)*SECDAY;
	    chi = puma_chi(puma, lnr, wr, rms);
	    chisrt[i+j*nrgrid] = chi;
	    chidx[i+j*nrgrid] = i+j*nrgrid;

	    if(sqrt(rms[2]*rms[2]+rms[3]*rms[3]) < chi0) {
	       RFIT = R;
	       VFIT = VR;
	       chi0 = sqrt(rms[2]*rms[2]+rms[3]*rms[3]);
	    }
	    if(PUMA_VERB > 1) {
	       printf("%3d %11.4e %7.3f  Chi: %9.3f  RMS: %9.3f %9.3f  Curv: %9.3f %9.3f %9.3f\n", 
		      i, R, VR, chi, rms[0], rms[1], rms[2], rms[3],
		      sqrt(rms[2]*rms[2]+rms[3]*rms[3]));
	    }

	    if( !(chi>0 || chi<=0)) chi = 1e38;	// fix any nans
	    if(JSON) {
	       printf("{ \"ir\": %2d, \"iv\": %2d, \"r_au\": %11.4e, \"logr\": %7.3f, \"v\": %7.3f, \"vr/r\": %7.3f, \"chi\": %12.4e}%s\n",
		      i, j, R, log(R)/log(10), VR, wr, chi,
		      (i==nrgrid-1&&j==nvgrid-1)?"],\n":",");

	    } else if(fpgrid != NULL) {
	       fprintf(fpgrid, "%3d %3d %11.4e %10.6f %7.3f %7.3f %12.4e\n",
		       i, j, R, log(R)/log(10), VR, wr, chi);

	    }
	 }
      }

/* Identify nrad acceptable R,VR, save as chirad and chivel */
      tsort_d(nrgrid*nvgrid, chisrt, chidx);
      for(k=nrad=0; k<nrgrid*nvgrid; k++) {
	 if(chisrt[k]-chisrt[0] < 1.0/(2*nfit-4)) {

	    i = chidx[k] % nvgrid;
	    j = chidx[k] / nvgrid;
	    R = rmin * exp((i*log(rmax/rmin))/(nrgrid-1));
	    VR = vmin + (j*(vmax-vmin))/(nvgrid-1);

	    chirad[nrad] = R * AU;
	    chivel[nrad++] = VR * 1e3;
	    if(PUMA_VERB > 1) {
	       printf("OK chi %.4f r= %.3e v= %7.3f   %d %d %d\n", 
		      chisrt[k], R, VR, 
		      chidx[k], chidx[k]/nvgrid, chidx[k]%nvgrid);
	    }
	 }
      }
      if(PUMA_VERB > 0) {
	 printf("Around R,VR %.3e %.3f are %d OK chi from %.4f to %.4f\n",
		RFIT,VFIT, nrad, chisrt[0], chisrt[nrad-1]);
      }

      lnr = log(RFIT);
      wr = VFIT*1e3/(RFIT*AU)*SECDAY;
   }
/////////////////////////////
///      GRID SEARCH      ///
/////////////////////////////



///////////////////////////////////
///      LEAST SQUARES FIT      ///
///////////////////////////////////
   
   if(dofit) {
      lnr = log(RFIT);
      wr = VFIT*1e3/(RFIT*AU)*SECDAY;
      chi0 = puma_fit(puma, &lnr, &wr, &dlnr, &dwr, &rwcov);

      RFIT = exp(lnr);			// [AU]
      VFIT = wr / SECDAY * RFIT;	// [km/s]

   }
///////////////////////////////////
///      LEAST SQUARES FIT      ///
///////////////////////////////////
   


/* Execute puma_chi to get GC fit for ephemeris points as well as fit points */
   chi0 = puma_chi(puma, lnr, wr, rms);

////////////////////////////////////////////////////////////////
   }		// for(iter=0; ...)
////////////////////////////////////////////////////////////////

/* Tell us all the results from the fits */
   if(PUMA_VERB > 0) puma_dump(puma, stdout);


/* Get the fit parameters */
   puma_fitpar(puma, &nfit, &ndet, &whichorigin, &ndate,
	       &t0, X, V, A, gcpole, gcorigin, gcphi);

/* Get residuals info for the reference observation */
   sigma = puma_resid(puma, whichorigin, &rafit, &decfit, &delta, &vrad, &wtan, &vpa, &xerr, &terr, &ragct, &decgct);

/* Radian and tangential velocity wrt observer [km/s] */
//   vrad0 = vrad;
//   vtan0 = (delta*AU*1e-3) * (wtan*atan(1)/45/SECDAY);


   if(JSON) {
/* Write the (ecliptic) state vector at the reference time */
      printf("\"refmjd\": [ %12.6f ], \"loc\": [ %9.6f, %9.6f, %9.6f ], \"vel\": [ %9.6f, %9.6f, %9.6f ],\n", 
	     t0, X[0]/AU, X[1]/AU, X[2]/AU, V[0]*1e-3, V[1]*1e-3, V[2]*1e-3);
   }




//////////////////////////////////
///      OUTPUT RESIDUALS      ///
//////////////////////////////////
   if(fpres != NULL || JSON) {
/* npt <= 0 -> no ephem MJDs, so emit fit resids only */
      if(JSON) {
	 printf("\"nobs\": %d, \"resid\": [\n", nfit);
      } else {
	 fprintf(fpres, "# Observations and residuals:\n");
	 fprintf(fpres, "# MJDobs      dti[day]  xresid    tresid      sigma     RAobs      Decobs     RAfit     Decfit   Distobs   wtan    PA     ID\n");
      }
      for(i=0; i<nfit; i++) {
	 /* Get residual info for this observation */
	 sigma = puma_resid(puma, i+1, &rafit, &decfit, &delta, &vrad, &wtan, &vpa, &xerr, &terr, &ragct, &decgct);
	 /* Get observation information from detection array */
	 puma_obsinfo(puma, i+1, &eph, &t0, &mjd, &raobs, &decobs);
	 puma_id(puma, i+1, 0, line);
	 if(JSON) {
	    printf("  { \"id\": \"%s\", \"MJDobs\": %11.6f, \"dt\": %9.4f, \"xresid\": %9.3f, \"tresid\": %9.3f, \"sigma\": %9.3f, \"RAobs\": %11.6f, \"Decobs\": %11.6f, \"RAfit\": %11.6f, \"Decfit\": %11.6f, \"Delta\": %8.5f, \"wtan\": %7.4f, \"wpa\": %7.2f}%s\n", 
		   line, mjd, mjd-t0, xerr, terr, sigma,
		   raobs, decobs, rafit, decfit, delta, wtan, vpa,
		   (i==nfit-1)?"],\n":",");
	 } else {
	    fprintf(fpres, "%12.6f %9.4f %9.3f %9.3f %9.3f %10.6f %10.6f %10.6f %10.6f %8.6f %8.4f %7.2f %s\n",
		    mjd, mjd-t0, xerr, terr, sigma,
		    raobs, decobs, rafit, decfit, delta, wtan, vpa, line);
	 }
      }
   }
//////////////////////////////////
///      OUTPUT RESIDUALS      ///
//////////////////////////////////





////////////////////////////////////////////
///      OUTPUT BEST-FIT EPHEMERIDES     ///
////////////////////////////////////////////
   if(fpeph != NULL) {
      if(JSON) {
	 printf("\"neph\": %d, \"ephem\": [\n", npt-nfit);
      } else {
	 fprintf(fpeph, "# MJDobs     dti[day]  x_GC-fit  t_GC-fit     theta     RAGCt      DecGCt     RAfit     Decfit   Distobs   wtan    PA     ID\n");
      }

      //      printf("\"neph\": %d, \"ephem\": [\n", npt-nfit);
      for(i=nfit; i<npt; i++) {
	 /* Get residual info for this ephemeris request */
	 sigma = puma_resid(puma, i+1, &rafit, &decfit, &delta, &vrad, &wtan, &vpa, &xerr, &terr, &ragct, &decgct);
	 /* Get time information from detection array */
	 puma_obsinfo(puma, i+1, &eph, &t0, &mjd, &raobs, &decobs);
	 puma_id(puma, i+1, 0, line);
	 if(JSON) {
	    printf("  { \"id\": \"%s\", \"MJDobs\": %11.6f, \"dt\": %9.4f, \"xresid\": %9.3f, \"tresid\": %9.3f, \"theta\": %9.3f, \"RAgct\": %11.6f, \"Decgct\": %11.6f, \"RAfit\": %11.6f, \"Decfit\": %11.6f, \"Delta\": %8.5f, \"wtan\": %7.4f, \"wpa\": %7.2f}%s\n", 
		   line, mjd, mjd-t0, xerr, terr, sigma,
		   ragct, decgct, rafit, decfit, delta, wtan, vpa,
		   (i==npt-1)?"],\n":",");
	 } else {

	    fprintf(fpeph, "%12.6f %9.4f %9.3f %9.3f %9.3f %10.6f %10.6f %10.6f %10.6f %8.6f %7.4f %7.2f %s\n",
		    mjd, mjd-t0, xerr, terr, sigma,
		    ragct, decgct, rafit, decfit, delta, wtan, vpa, line);
	 }
      }
   }
////////////////////////////////////////////
///      OUTPUT BEST-FIT EPHEMERIDES     ///
////////////////////////////////////////////



//////////////////////////////////
///      CLOSEST APPROACH      ///
//////////////////////////////////
   rperi = puma_perigee(puma, &mjdperi, &raperi, &decperi, &wperi);
//////////////////////////////////
///      CLOSEST APPROACH      ///
//////////////////////////////////








////////////////////////////////////////////
///      OUTPUT STATISTICAL SCATTER      ///
////////////////////////////////////////////
/* Compute a statistical scatter for each ephemeris MJD? */
   if(fpscat != NULL) {

      if(JSON) {
//	 printf("\"nscat\": %d, \"scatter\": [\n", NephMC);
	 printf("\"nscat\": %d, \"scatter\": [", NephMC);
      } else {
	 fprintf(fpscat, "# MC  eph    MJD           RA        Dec      chi/N     R[AU]  v[km/s]\n");
      }

/* Sample points that have consistent R,VR and unit vectors */
      for(i=k=0; i<NephMC; i++) {
	 slnr = lnr;
	 swr = wr;
	 chi = puma_perturb(puma, &slnr, &swr, dlnr, dwr, rms);
	 R = AU * exp(slnr);
	 VR = swr * R / SECDAY;

/* If the solution is sufficiently OK write it */
	 if((chi-chi0) < chimax) {
	    for(j=nfit; j<npt; j++) {
/* Return information from detection array */
	       puma_obsinfo(puma, j+1, &eph, &t0, &mjd, &raobs, &decobs);

/* Get the residuals */
	       sigma = puma_resid(puma, j+1, &rafit, &decfit, &delta, &vrad, &vpa, &wtan, &xerr, &terr, &ragct, &decgct);
	       if(JSON) {
		  printf("%s\n  { \"mc\": %3d,  \"eph\": %2d, \"MJDobs\": %12.6f, \"RA\": %11.6f, \"Dec\": %11.6f, \"chi\": %7.3f, \"r_au\": %8.5f, \"vr\": %7.4f}", 
			 (k==0)?"":",",
			 i, j+1, mjd, rafit, decfit, chi, R/AU, VR/1e3);
		  k++;
	       } else {
		  fprintf(fpscat, "%4d %3d %12.6f %9.5f %9.5f %11.3e %8.6f %7.2f\n", 
			  i, j+1, mjd, rafit, decfit, chi, R/AU, VR/1e3);
	       }
	    }
	 }
      }
      if(JSON) {
	 printf("],\n");
      }
   }
////////////////////////////////////////////
///      OUTPUT STATISTICAL SCATTER      ///
////////////////////////////////////////////



////////////////////////////////////////
///      OUTPUT SUMMARY RESULTS      ///
////////////////////////////////////////

/* Return selected Keplerian orbit parameters wrt Sun or Earth */

/* Test whether Earth-bound */
   geobound = 1;
   puma_kepler(puma, geobound, &Egeo, &ecc, &incl, &a, &period);
   E = Egeo;

/* Nope, not bound to Earth */
   if(E > 0) {
      geobound = 0;
      puma_kepler(puma, geobound, &Esun, &ecc, &incl, &a, &period);
      E = Esun;
   }

/* Total velocity wrt Sun */
//   vsun = sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]);

/* Uncertainty in energy with respect to Earth [ (km/s)^2 ] */
   dEgeo = geo_E_error(lnr+log(AU) /*ln(m)*/, dlnr,
                       wr/SECDAY /*sec^-1*/, dwr/SECDAY /*sec^-1*/,
                       rwcov /*corcoeff*/, Egeo*1e6 /*(m/s)^2*/) * 1e-6;

   if(JSON) {
      printf("\"results\": {\n");
      printf("  \"bound-to\": \"%s\",\n", geobound==0?"sun":"earth");
      printf("  \"E\": %.3f,\n", E);
      printf("  \"P\": %.3f,\n", period);
      printf("  \"a\": %.3f,\n", a);
      printf("  \"ecc\": %.3f,\n", ecc);
      printf("  \"i\": %.3f,\n", incl);
      printf("  \"Egeo\": %.3f,\n", Egeo);
      printf("  \"dEgeo\": %.3f,\n", dEgeo);

   } else {
// JTEST
//      double cc[5], c = puma_const(puma, &cc);
//      printf("chin const %.3f %.5f\n", c, cc[0]);
// JTEST
      printf("Orbit with respect to the %s\n", geobound==0?"Sun":"Earth");
      printf("E,P,a,ecc,i,Egeo,dEgeo= %7.2f %6.3f %6.3f %5.3f %6.3f %6.3f %6.3f",
             E, period, a, ecc, incl, Egeo, dEgeo);
      printf(" // [(km/s)^2] [%s] [%s] [] [deg]\n",
	     geobound==0?"year":"day", geobound==0?"AU":"km");
   }

   if(JSON) {
      printf("  \"ndet\": %d,\n", nfit);
      printf("  \"ndate\": %d,\n", ndate);
      printf("  \"t0_mjd\": %.6f,\n", t0);
      printf("  \"r_km\": %.0f,\n",  1e-3*AU*exp(lnr));
      printf("  \"r_au\": %.5f,\n", exp(lnr));
      printf("  \"dlnr\": %.4f,\n", dlnr);
      printf("  \"Vr\": %.5f,\n", wr/SECDAY*exp(lnr)*AU*1e-3);
      printf("  \"Vr/r\": %.5f,\n", wr);
      printf("  \"dVr\": %.5f,\n", dwr);
//      printf("  \"Vsun\": %.5f,\n", vsun*1e-3);
//      printf("  \"Vrad\": %.5f,\n", vrad0);
//      printf("  \"Vtan\": %.5f,\n", vtan0);
      printf("  \"rvcov\": %.3f,\n", rwcov);
      printf("  \"chi/N\": %.3f,\n", chi0);
      printf("  \"xrms\": %.3f,\n", rms[0]);
      printf("  \"trms\": %.3f,\n", rms[1]);
      printf("  \"xcrv\": %.3f,\n", rms[2]);
      printf("  \"tcrv\": %.3f,\n", rms[3]);
      printf("  \"tperi_mjd\": %.2f,\n", mjdperi);
      printf("  \"rperi_au\": %.6f,\n", rperi);
      printf("  \"wperi\": %.4f,\n", wperi);
      printf("  \"raperi\": %.1f,\n", raperi);
      printf("  \"decperi\": %.1f,\n", decperi);
      printf("  \"units\": { \"E\": \"[(km/s)^2]\",  \"Egeo\": \"[(km/s)^2]\", \"P\": \"[%s]\", \"a\": \"[%s]\", \"i\": \"[deg]\" }\n",
	     geobound==0?"year":"day", geobound==0?"AU":"km");

      printf("}\n}\n");

     
   } else {
      printf("ndet=  %12d  // number of detections fitted\n", nfit);
      printf("ndate= %12d  // number of distinct epochs\n", ndate);
      printf("t0=    %12.6f  // [MJD] reference time (observer)\n", t0);
      printf("r=     %12.0f  // [km] distance from barycenter to obj at t=ref\n", 1e-3*AU*exp(lnr));
      printf("rau=   %12.4e  // [AU] distance from barycenter to obj at t=ref\n", exp(lnr));
      printf("dlnr=  %12.4f  // uncertainty in ln(r)\n", dlnr);
      printf("Vr=    %12.5f  // [km/s] radial velocity of obj wrt barycenter\n", wr/SECDAY*exp(lnr)*AU*1e-3);
      printf("Vr/r=  %12.5f  // [/day] radial vel/dist of obj wrt barycenter\n", wr);
      printf("dV/r=  %12.5f  // [/day] uncertainty in Vr/r\n", dwr);
//      printf("Vsun=  %12.5f  // [km/s] total velocity of obj wrt sun\n", vsun*1e-3);
//      printf("Vrad=  %12.5f  // [km/s] radial velocity of obj wrt observer\n", vrad0);
//      printf("Vtan=  %12.5f  // [km/s] tangential velocity of obj wrt observer\n", vtan0);
      printf("rvcov= %12.3f  // covariance between ln(r) and vr/r\n", rwcov);
      printf("chi/N= %12.3f  // chi^2/Ndof\n", chi0);
      printf("xrms=  %12.3f  // [arcsec] cross-track RMS\n", rms[0]);
      printf("trms=  %12.3f  // [arcsec] tangential RMS\n", rms[1]);
      printf("xcrv=  %12.3e  // [arcsec/day^2] cross-track curvature\n", rms[2]);
      printf("tcrv=  %12.3e  // [arcsec/day^2] tangential curvature\n", rms[3]);
//   printf("xcrvt2=%12.3f  // [arcsec] cross-track curv times <dt>^2\n", rms[2]*dt*dt);
//   printf("tcrvt2=%12.3f  // [arcsec] tangential curv times <dt>^2\n", rms[3]*dt*dt);
      printf("tperi= %12.2f  // [MJD] time of closest approach\n", mjdperi);
      printf("rperi= %12.6f  // [AU] distance of closest approach\n", rperi);
      printf("wperi= %12.4f  // [deg/day] ang vel at closest approach\n", wperi);
      printf("raperi=  %10.1f  // [deg] RA at closest approach\n", raperi);
      printf("decperi= %10.1f  // [deg] Dec at closest approach\n", decperi);
   }
////////////////////////////////////////
///      OUTPUT SUMMARY RESULTS      ///
////////////////////////////////////////




/* Release memory */
   puma_free(puma);

   if(niter > 0 || PUMA_VERB > 0) {
/* End time */
      gettimeofday(&tv1, NULL);
      dt = (tv1.tv_sec - tv0.tv_sec) + 1e-6*(tv1.tv_usec - tv0.tv_usec);

      printf("Time= %.6f  Nxrock= %d\n", dt, NXROCK);
   }

   exit(0);
}





// Input options:
//
// MJD RA Dec (poserr tngerr lng lat elev from command line, obsID='Anon')
// MJD RA Dec poserr (tngerr=poserr lng lat elev from command line, obsID='Anon')
// MJD RA Dec poserr tngerr (lng lat elev from command line)
// MJD RA Dec poserr tngerr lng lat elev  (obsID='Anon')
// MJD RA Dec poserr tngerr lng lat elev obsID
//
// Units:
//   MJD [MJD]
//   RA,Dec [deg]
//   poserr,tngerr [arcsec] = [~FWHM/SNR]  (Use poserr<0 for eval only)
//   lng,lat [deg]
//   elev [m]



/* Read an input file, return number of detections and data */
int rock_input(char *infile,
               double lng_default, double lat_default, double elev_default,
               double lng_eph, double lat_eph, double elev_eph,
               double poserr_default, char *epharg, void *p)
{
   int i, n, nitem, neph;
   double mjd, ra, dec, poserr, tngerr, lng, lat, elev;
   char line[1024], obsid[1024];
   FILE *fp;

/* Open input */
   if(strcmp(infile, "-") == 0) {
      fp = stdin;
   } else if( (fp = fopen(infile, "r")) == NULL) {
      fprintf(stderr, "Cannot open %s for reading\n", infile);
      exit(1);
   }

/* Read all the data */
   for(n=0;   ; n++) {
      if(fgets(line, 1024, fp) == NULL) break;
      if(line[0] == '#') {n--; continue;}	/* Skip comments */
      nitem = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %s",
		     &mjd, &ra, &dec, &poserr, &tngerr, &lng, &lat, &elev,
		     obsid);
      if(nitem <= 0) {n--; continue;}
      if(nitem < 3) {
	 fprintf(stderr, "Error reading from `%s' at line %d\n", infile, n);
	 fprintf(stderr, "Line `%s'\n", line);
	 exit(1);
      }
      if(nitem < 9) strcpy(obsid, "Anon");
      if(nitem < 8) {
	 if(lng_default < -999 || lat_default < -999 || elev_default < -999) {
	    fprintf(stderr, "Must specify observer -lng -lat -elev\n");
	    exit(1);
	 }
	 lng = lng_default;
	 lat = lat_default;
	 elev = elev_default;
      }
      if(nitem < 4) poserr = poserr_default;
      if(nitem < 5) tngerr = poserr;

      i = puma_obs(p, 0, mjd, ra, dec, poserr, tngerr, lng, lat, elev);

/* Patch up the id's */
      if(strcmp(obsid, "Anon") == 0) sprintf(obsid, "%d", i);
      puma_id(p, i, 1, obsid);
   }
   if(strcmp(infile, "-") != 0) fclose(fp);

   double *ephmjd;
   double *ephlng;
   double *ephlat;
   double *ephelev;
   char **ephid;

   ephmjd = (double *)calloc(MAXEPHMJD, sizeof(double));
   ephlng = (double *)calloc(MAXEPHMJD, sizeof(double));
   ephlat = (double *)calloc(MAXEPHMJD, sizeof(double));
   ephelev = (double *)calloc(MAXEPHMJD, sizeof(double));
   ephid = (char **)calloc(MAXEPHMJD, sizeof(char *));

   for(i=0; i<MAXEPHMJD; i++) ephlng[i] = ephlat[i] = ephelev[i] = -1000;

/* Extra MJD(s) for evaluation? */
   neph = parse_eph(epharg, MAXEPHMJD, ephmjd, ephid, ephlng, ephlat, ephelev);

   for(i=0; i<neph; i++) {
/* No site data from the file provided? */
      if(ephlng[i] < -999 || ephlat[i] < -999 || ephelev[i] < -999) {
/* Site data from command line? */
         if(lng_eph > -999 && lat_eph > -999 && elev_eph > -999) {
            ephlng[i] = lng_eph;
            ephlat[i] = lat_eph;
            ephelev[i] = elev_eph;
/* Just use the site info from the last observation */
         } else {
            ephlng[i] = lng;
            ephlat[i] = lat;
            ephelev[i] = elev;
         }
      }

      n = puma_obs(p, 1, (double)ephmjd[i], 0., 0., 0., 0.,
                   ephlng[i], ephlat[i], ephelev[i]);
      puma_id(p, n, 1, ephid[i]);
   }
   
   return(n);
}

/* Parse a CSV string with real variables, return array and number */
int parse_eph(char *arg, int maxvar, double *mjd, char **ephid,
              double *lng, double *lat, double *elev)
{
   int n, narg;
   char *comma;
   double t0, t1, dt;
   FILE *fp;
   #define LINESZ 1024
   char line[LINESZ];

   if(arg == NULL) return(0);
/* Is this "file:filename" containing a list of MJDs? */
   if (strncmp(arg, "file:", 5) == 0) {
      if (strlen(arg) < 6) {
         fprintf(stderr, "weird ephmjd file %s\n", arg);
         exit(1);
      }
      if((fp = fopen(&arg[5], "r")) == NULL) {
         fprintf(stderr, "Cannot open ephmjd file %s for reading\n", &arg[5]);
         exit(1);
      }
      n = 0;
      while (fgets(line, LINESZ, fp) > 0) {
         ephid[n] = malloc(32);
         if ( (narg=sscanf(line, "%lf %s %lf %lf %lf", &mjd[n], ephid[n],
                           &lng[n], &lat[n], &elev[n])) < 1) {
            break;
         }
	 if(narg == 1) snprintf(ephid[n], 32, "%d", n+1);
         n++;
         if (n == MAXEPHMJD) {
            fprintf(stderr, "max %d ephmjds reached\n", MAXEPHMJD);
            break;
         }
      }
      fclose(fp);

/* Look for a colon: if so parse as start:end:interval */
   } else if(index(arg, ':') != NULL) {
      if(sscanf(arg, "%lf:%lf:%lf", &t0, &t1, &dt) < 3) {
	 fprintf(stderr, "Cannot parse t0:t1:dt from %s\n", arg);
	 return(0);
      }
      if(strcmp(arg+strlen(arg)-1, "m") == 0) dt /= 1440;
      if(strcmp(arg+strlen(arg)-1, "h") == 0) dt /= 24;
      for(n=0; n<=MIN(maxvar-1, (t1-t0)/dt); n++) {
	 mjd[n] = t0 + n*dt;
         ephid[n] = malloc(32);
	 snprintf(ephid[n], 32, "%d", n+1);
      }

/* Parse as a vanilla CSV */
   } else {
      for(n=0; n<maxvar; n++) {
	 if(sscanf(arg, "%lf", &mjd[n]) <= 0) break;
         ephid[n] = malloc(32);
	 snprintf(ephid[n], 32, "%d", n+1);
	 if( (comma = index(arg, ',')) == NULL) {n++; break;}
	 arg = comma + 1;
	 if(strlen(arg) == 0) {n++; break;}
      }
   }
   return(n);
}

#define GMEARTH (3.986004418e14) /* [m^3/s^2] */

/* Uncertainty in energy with respect to Earth */
/* x = lnr = log(r/m), y = wr = vr/r [/sec], */
/* E = exp(2*x)*(y^2+omega^2)/2 - exp(-x)*GM [(m/s)^2]*/
double geo_E_error(double lnr, double dlnr, double wr, double dwr,
                   double rwcov, double E)
{
   double dE, r=exp(lnr), det, t2, alpha, c, s, dEdx, dEdy, du2, dv2;

/* Limit stupid values */
   if(dlnr < 1e-30) dlnr = 1e-30;
   if(dwr < 1e-30) dwr = 1e-30;

/* Determinant of covariance matrix -> inverse covariance matrix */
   det = dlnr*dlnr * dwr*dwr * (1 - rwcov*rwcov);

/* Rotation angle to diagonalize inverse covariance matrix from x,y to u,v */
   t2 = 1e20;
   if(dlnr != dwr) t2 = 2*dlnr*dwr*rwcov / (dwr*dwr-dlnr*dlnr);
   alpha = 0.5 * atan(t2);
   c = cos(alpha);
   s = sin(alpha);

/* Derivative of E with respect to x,y */
   dEdx = 2*E + 3*GMEARTH/r;
   dEdy = r*r * wr;

/* Uncertainties in rotated, diagonal, u,v coordinates */
   du2 = det / (c*c*dwr*dwr + 2*s*c*dlnr*dwr*rwcov + s*s*dlnr*dlnr);
   dv2 = det / (c*c*dlnr*dlnr - 2*s*c*dlnr*dwr*rwcov + s*s*dwr*dwr);

/* (dE/dx dx/du + dE/dy dy/du)^2 * du^2 + (dE/dx dx/dv + dE/dy dy/dv)^2 * dv^2 */
   dE = (dEdx*c+dEdy*s)*(dEdx*c+dEdy*s) * du2 +
       (-dEdx*s+dEdy*c)*(dEdx*c+dEdy*s) * dv2;
   if(dE > 0) dE = sqrt(dE);    /* should always be positive! */

   return(dE);
}
