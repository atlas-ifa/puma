/* Convert MPC rho*cos(phi') rho*sin(phi') to lat and height above ellipsoid */
/* 20200107 - John Tonry implementation of Lin & Wang 1995 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* The WGS84 model of the Earth */
#define WGS84_a (6378137.0)		/* [m] */
#define WGS84_finv (298.257223563)	/* [] inverse flattening */

double a=WGS84_a;
double b=WGS84_a*(1-1/WGS84_finv);
double f=1/WGS84_finv;


double fm(double p, double Z, double m)
{
   return(p*p/(a+2*m/a)/(a+2*m/a) +
	  Z*Z/(b+2*m/b)/(b+2*m/b) - 1);
}

double dfm(double p, double Z, double m)
{
   return(-4*(p*p/a/(a+2*m/a)/(a+2*m/a)/(a+2*m/a) +
	      Z*Z/b/(b+2*m/b)/(b+2*m/b)/(b+2*m/b)));
}



int main(int argc, char **argv)
{
   int tompc=0;
   double rhocos, rhosin, lat, h;
   double m, p, Z, pe, Ze, dr=atan(1)/45;

   if(argc < 3) {
      printf("Syntax: mpclat rhocos rhosin -> lat[deg] elev[m]\n");
      printf("  converts MPC geocentric coords to geodetic latitude\n");
      printf("e.g. mpclat  0.936235 +0.351547 > 20.707571 3040.9\n");
      printf("\nThe inverse function also exists:\n");
      printf("  mpclat lat[deg] elev[m] tompc -> rhocos rhosin\n");
      printf("e.g. mpclat 20.707571 3040.9 tompc > 0.936235  0.351547\n");
      exit(1);
   }

   sscanf(argv[1], "%lf", &rhocos);
   sscanf(argv[2], "%lf", &rhosin);

/* Inverse function to convert lat,elev to MPC cosphi sinphi */
   if(argc>3 && strcmp(argv[3], "tompc") == 0) tompc = 1;

   if(tompc == 0) {
      p = a * rhocos;
      Z = a * rhosin;
      m = (a*b*pow(a*a*Z*Z+b*b*p*p,1.5)-a*a*b*b*(a*a*Z*Z+b*b*p*p)) /
	 (2*(a*a*a*a*Z*Z+b*b*b*b*p*p));
      m -= fm(p,Z,m) / dfm(p,Z,m);
      pe = p / (1+2*m/a/a);
      Ze = Z / (1+2*m/b/b);

      lat = atan(a*a*Ze/b/b/pe) / dr;
      h = sqrt((p-pe)*(p-pe)+(Z-Ze)*(Z-Ze));
      printf("%.6f %.1f\n", lat, h);
      return(0);
   }

/* Convert lat,elev to MPC cosphi,sinphi */
   double C, S, rho, phi;

   lat = rhocos * dr;
   h = rhosin;

   C = 1/sqrt(cos(lat)*cos(lat)+(1-f)*(1-f)*sin(lat)*sin(lat));
   S = (1-f)*(1-f) * C;
   p = (a*C+h) * cos(lat);
   Z = (a*S+h) * sin(lat);
   rho = sqrt(p*p+Z*Z);
   phi = atan2(Z, p);
   rhocos = rho/a * cos(phi);
   rhosin = rho/a * sin(phi);
   printf("%9.6f %9.6f\n", rhocos, rhosin);
   return(0);

}
