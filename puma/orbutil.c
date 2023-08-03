/* Utilities to integrate an orbit */
/* v1.0 - 151127 John Tonry, broken out of orbit.c */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "orbit.h"

extern int VERB;

static int DTVERB=0;	/* Gabby messages from tripdt? */

void lp_xyzmoon(REAL mjd, VEC *moon);
void lp_xyzsun(REAL mjd, VEC *moon);

/* Local sidereal time (rad) at modified Julian date mjd at E long (rad) */
/* From Astronomy on the Personal Computer */
// FIXME?  UT1 vs MJD(UTC), argh
REAL lst(REAL mjd, REAL elng)
{
   REAL mjd0, ut, T0, T, gmst, lst;
   REAL pi=4*ATAN(1.0);
   mjd0 = floor(mjd);
   ut = SECDAY * (mjd-mjd0);		// [sec]
   T0 = (mjd0 - 51544.5) / 36525.0;	// [century]
   T = (mjd - 51544.5) / 36525.0;	// [century]
   gmst = 24110.54841 + 8640184.812866 * T0 + 1.0027379093 * ut
      + (0.093104 - 6.2e-6*T)*T*T;	// [sec]
   gmst *= 2*pi/SECDAY;		// [rad]
//   printf("GMST= %8.1f %8.3Lf\n", ut, fmod(gmst,2*pi)*180/pi);
   lst = gmst + elng;			// [rad]
   return(fmod(lst, 2*pi));		// [Rad]
}

// http://www.gpstk.org/doxygen/classgpstk_1_1GeodeticFrames.html

// Reference: IERS Technical Note 21, IERS Conventions (1996), 
// Dennis D. McCarthy, U.S. Naval Observatory.
// 
// The conventional terrestrial system (CTS) or Earth-centered,
// Earth-fixed frame (ECEF), is related to the conventional inertial
// system (CIS) by four things: a) precession and b) nutation of the
// Earth and c) Earth rotation and d) polar motion. The transformation
// between a vector X(CTS) in the terrestrial (ECEF) frame and the vector
// X(CIS) in the inertial frame is
// 
//  X(CIS) = P * N * R * W * X(CTS)
//     where
//  W is the transformation using Earth Orientation Parameters
//     xp, yp (pole coordinates obtained from the IERS bulletin),
//  R is the effect of both Earth rotation and
//     precession and nutation in the right ascension,
//  N is the nutation matrix,
//  P is the precession matrix.
//  Reference: IERS Technical Note 21, IERS Conventions (1996), Chapter 5.
//  IF
//    R1(a) =  [ 1    0      0    ]
//             [ 0  cos(a) sin(a) ]
//             [ 0 -sin(a) cos(a) ]
//    R2(b) =  [ cos(b) 0 -sin(b) ]
//             [   0    1    0    ]
//             [ sin(b) 0  cos(b) ]
//    R3(c) =  [  cos(c) sin(c) 0 ]
//             [ -sin(c) cos(c) 0 ]
//             [    0      0    1 ]
//  and if
//    T = (t-t0)(in days)/36525.0 days
//    where
//    t0 = J2000 = January 1 2000 12h UT = 2451545.0JD
//  THEN ----------------------------------------------------------------
//   [PRECESSION IAU76]
//    P = R3(zeta)*R2(-theta)*R3(z)
//   where
//    zeta  = 2306.2181*T + 0.30188*T^2 + 0.017998*T^3 seconds of arc
//    theta = 2004.3109*T - 0.42665*T^2 - 0.041833*T^3 seconds of arc
//    z     = 2306.2181*T + 1.09468*T^2 + 0.018203*T^3 seconds of arc
//  AND -----------------------------------------------------------------
//   [NUTATION IAU76]
//    N = R1(-eps)*R3(dpsi)*R1(eps+deps)
//   where
//    eps  = obliquity of the ecliptic
//    deps = nutation in obliquity
//    dpsi = nutation in longitude (counted in the ecliptic)
//    eps = (84381.448 - 46.8150*T - 0.00059*T^2 +0.001813*T^3) seconds of arc
//   and the principal terms in the series for the other two parameters are
//   ( " denotes seconds of arc)
//   (IAU76)
//    deps   =
//      ( 9.205356 + 0.000886*T)*cos(Omega)"       + 0.001553*sin(Omega)"
//    + ( 0.573058 - 0.000306*T)*cos(2F-2D+2Omega)"- 0.000464*sin(2F-2D+2Omega)"
//    + ( 0.097864 - 0.000048*T)*cos(2F+2Omega)"   + 0.000136*sin(2F+2Omega)"
//    + (-0.089747 + 0.000047*T)*cos(2Omega)"      - 0.000029*sin(2Omega)"
//    + ( 0.007388 - 0.000019*T)*cos(-L')"         + 0.000198*sin(-L')"
//    + ( 0.022440 - 0.000068*T)*cos(Larg)"        - 0.000018*sin(Larg)"
//    + (-0.000687 + 0.000000*T)*cos(L)"           - 0.000039*sin(L)"
//    etc...
//    dpsi   =     
//      (-17.206277- 0.017419*T)*sin(Omega)"       + 0.003645*cos(Omega)"
//    + (-1.317014 - 0.000156*T)*sin(2F-2D+2Omega)"- 0.001400*cos(2F-2D+2Omega)"
//    + (-0.227720 - 0.000023*T)*sin(2F+2Omega)"   + 0.000269*cos(2F+2Omega)"
//    + ( 0.207429 + 0.000021*T)*sin(2Omega)"      - 0.000071*cos(2Omega)"
//    + (-0.147538 + 0.000364*T)*sin(-L')"         + 0.001121*cos(-L')"
//    + (-0.051687 + 0.000123*T)*sin(Larg)"        - 0.000054*cos(Larg)"
//    + ( 0.071118 + 0.000007*T)*sin(L)"           - 0.000094*cos(L)"
//    etc...
//  OR
//   (IERS 1980)
//    deps   =   ( 9.2025 + 0.00089*T)*cos(Omega)"
//             + ( 0.5736 - 0.00031*T)*cos(2F-2D+2Omega)"
//             + ( 0.0977 - 0.00005*T)*cos(2F+2Omega)"
//             + (-0.0895 + 0.00005*T)*cos(2Omega)"
//             + ( 0.0054 - 0.00001*T)*cos(-L')"
//             + (-0.0007 + 0.00000*T)*cos(L)"
//    etc...
//    dpsi   =   (-17.1996- 0.01742*T)*sin(Omega)"
//             + (-1.3187 - 0.00016*T)*sin(2F-2D+2Omega)"
//             + (-0.2274 - 0.00002*T)*sin(2F+2Omega)"
//             + ( 0.2062 + 0.00002*T)*sin(2Omega)"
//             + (-0.1426 + 0.00034*T)*sin(-L')"
//             + ( 0.0712 + 0.00001*T)*sin(L)"
//    etc...
//   with
//    Larg    = L'+2F-2D+2Omega
//    Omega   = mean longitude of the lunar ascending node
//            = 125.04455501 degrees - 6962890.2665"*T + 7.4722"*T^2
//                + 0.007702"*T^3 - 0.00005939"*T^4
//    D       = mean elongation of the moon from the sun
//            = 297.85019547 degrees + 1602961601.2090"*T - 6.3706"*T^2
//                + 0.006593"*T^3 - 0.00003169"*T^4
//    F       = mean longitude of the moon - Omega
//            = 93.27209062 degrees + 1739527262.8478"*T - 12.7512"*T^2
//                + 0.001037"*T^3 + 0.00000417"*T^4
// 
//    L'      = mean anomaly of the sun
//            = 357.52910918 degrees + 129596581.0481"*T - 0.5532"*T^2
//                + 0.000136"*T^3 - 0.00001149"*T^4
//    L       = mean anomaly of the moon
//            = 134.96340251 degrees + 1717915923.2178"*T + 31.8792"*T^2
//                +0.051635"*T^3 - 0.00024470"*T^4
//  AND -----------------------------------------------------------------
//    R = R3(-GAST)
//    GAST = Greenwich hour angle of the true vernal equinox
//    GAST = Greenwich Apparent Sidereal Time
//    GAST = GMST + dpsi*cos(eps) + 0.00264"*sin(Omega) +0.000063"*sin(2*Omega)
//       (these terms account for the accumulated precession and nutation in
//          right ascension and minimize any discontinuity in UT1)
//    GMST = Greenwich hour angle of the mean vernal equinox
//         = Greenwich Mean Sidereal Time
//         = GMST0 + r*[(UT1-UTC)+UTC]
//    r    = is the ratio of universal to sidereal time
//         = 1.002737909350795 + 5.9006E-11*T' - 5.9e-15*T'^2
//    T'   = days'/36525
//    days'= number of days elapsed since t0 = +/-(integer+0.5)
//       and
//    (UT1-UTC) is taken from the IERS bulletin (seconds)
//    GMST0 = GMST at 0h UT1
//         = 6h 41m (50.54841+8640184.812866*T'+0.093104*T'^2-6.2E-6*T'^3)s
//  AND -----------------------------------------------------------------
//    W = R1(yp)*R2(xp)
//    where xp and yp are the polar coordinates of the celestial ephemeris pole,
//       taken from the IERS bulletin. (NB in the bulletin they are in units
//       of arcseconds, and they must be converted to radians by multiplying
//       by pi/180/3600.)


/* Rotate by phi [rad] around x axis: coords go CCW, vector goes CW */
void rotx(REAL phi, REAL *R1)
{
   R1[0 + 0*3] = 1;
   R1[1 + 0*3] = 0;
   R1[2 + 0*3] = 0;
   R1[0 + 1*3] = 0;
   R1[1 + 1*3] = COS(phi);
   R1[2 + 1*3] = SIN(phi);
   R1[0 + 2*3] = 0;
   R1[1 + 2*3] = -SIN(phi);
   R1[2 + 2*3] = COS(phi);
}

/* Rotate by phi [rad] around y axis: coords go CCW, vector goes CW */
void roty(REAL phi, REAL *R2)
{
   R2[0 + 0*3] = COS(phi);
   R2[1 + 0*3] = 0;
   R2[2 + 0*3] = -SIN(phi);
   R2[0 + 1*3] = 0;
   R2[1 + 1*3] = 1;
   R2[2 + 1*3] = 0;
   R2[0 + 2*3] = SIN(phi);
   R2[1 + 2*3] = 0;
   R2[2 + 2*3] = COS(phi);
}

/* Rotate by phi [rad] around z axis: coords go CCW, vector goes CW */
void rotz(REAL phi, REAL *R3)
{
   R3[0 + 0*3] = COS(phi);
   R3[1 + 0*3] = SIN(phi);
   R3[2 + 0*3] = 0;
   R3[0 + 1*3] = -SIN(phi);
   R3[1 + 1*3] = COS(phi);
   R3[2 + 1*3] = 0;
   R3[0 + 2*3] = 0;
   R3[1 + 2*3] = 0;
   R3[2 + 2*3] = 1;
}

/* Multiply matrices: C = A * B */
void matmult(REAL *A, REAL *B, REAL *C)
{
   REAL TMP[9];
   TMP[0] = A[0]*B[0] + A[1]*B[3] + A[2]*B[6];
   TMP[1] = A[0]*B[1] + A[1]*B[4] + A[2]*B[7];
   TMP[2] = A[0]*B[2] + A[1]*B[5] + A[2]*B[8];
   TMP[3] = A[3]*B[0] + A[4]*B[3] + A[5]*B[6];
   TMP[4] = A[3]*B[1] + A[4]*B[4] + A[5]*B[7];
   TMP[5] = A[3]*B[2] + A[4]*B[5] + A[5]*B[8];
   TMP[6] = A[6]*B[0] + A[7]*B[3] + A[8]*B[6];
   TMP[7] = A[6]*B[1] + A[7]*B[4] + A[8]*B[7];
   TMP[8] = A[6]*B[2] + A[7]*B[5] + A[8]*B[8];
   memcpy(C, TMP, 9*sizeof(REAL));
}

//  X(CIS=Inertial) = PNR * X(CTS=Terrestrial)
/* Return precession, nutation, and rotation matrix from mjd to 2000.0 */
void precnutrot(REAL mjd, REAL dUT1, REAL *PNR)
{
   REAL pi=4*ATAN(1.0), secrad=ATAN(1.0)/(45*3600);
   REAL T, zeta, theta, z;
   REAL P[9], N[9], R[9], TMP[9];
   REAL Omega, D, F, Lp, L;
   REAL eps, deps, dpsi;
   REAL Tp, ut1, r, GMST0, GMST, GAST;

   T = (mjd-51544.5) / 36525.0;

/* Precession */
   zeta  = T * (2306.2181 + T*(0.30188 + T*0.017998)) * secrad;
   theta = T * (2004.3109 + T*(-0.42665 - T*0.041833)) * secrad;
   z     = T * (2306.2181 + T*(1.09468 + T*0.018203)) * secrad;
   rotz(zeta, P);
   roty(-theta, TMP);
   matmult(P, TMP, P);
   rotz(z, TMP);
   matmult(P, TMP, P);

#ifdef DEBUG1
   printf("%s  %8.3Lf %8.3Lf %8.3Lf\n", "a", zeta*180/pi, theta*180/pi, z*180/pi);
   printf("%s  %9.5Lf %9.5Lf %9.5Lf\n", "P", P[0], P[1], P[2]);
   printf("%s  %9.5Lf %9.5Lf %9.5Lf\n", "P", P[3], P[4], P[5]);
   printf("%s  %9.5Lf %9.5Lf %9.5Lf\n", "P", P[6], P[7], P[8]);
#endif

/* Nutation */
   Omega = (3600*125.04455501 + T*(-6962890.2665 + 
            T*(7.4722 + T*(0.007702 - T*0.00005939)))) * secrad;
   D = (3600*297.85019547 + T*(1602961601.2090 + 
            T*(-6.3706 + T*(0.006593 - T*0.00003169)))) * secrad;
   F = (3600*93.27209062 + T*(1739527262.8478 + 
            T*(-12.7512 + T*(0.001037 + T*0.00000417)))) * secrad;
   Lp = (3600*357.52910918 + T*(129596581.0481 +
	    T*(-0.5532 + T*(0.000136 + T*-0.00001149)))) * secrad;
   L = (3600*134.96340251 + T*(1717915923.2178 + 
            T*(31.8792 + T*(0.051635 - T*0.00024470)))) * secrad;

   eps = (84381.448 + T*(-46.8150 + T*(-0.00059 + T*0.001813))) * secrad;

   deps = ( (9.2025 + 0.00089*T) * COS(Omega) + 
	    (0.5736 - 0.00031*T) * COS(2*F-2*D+2*Omega) + 
	    (0.0977 - 0.00005*T) * COS(2*F+2*Omega) +
	    (-0.0895 + 0.00005*T) * COS(2*Omega) +
	    (0.0054 - 0.00001*T) * COS(-Lp) +
	    (-0.0007 + 0.00000*T) * COS(L) ) * secrad;

   dpsi = ( (-17.1996 - 0.01742*T) * SIN(Omega) +
	    (-1.3187 - 0.00016*T) * SIN(2*F-2*D+2*Omega) +
	    (-0.2274 - 0.00002*T) * SIN(2*F+2*Omega) +
	    ( 0.2062 + 0.00002*T) * SIN(2*Omega) +
	    (-0.1426 + 0.00034*T) * SIN(-Lp) +
	    ( 0.0712 + 0.00001*T) * SIN(L) ) * secrad;

   rotx(-eps, N);
   rotz(dpsi, TMP);
   matmult(N, TMP, N); 
   rotx(eps+deps, TMP);
   matmult(N, TMP, N);

#ifdef DEBUG1
   printf("%s  %9.5Lf %9.5Lf %9.5Lf\n", "N", N[0], N[1], N[2]);
   printf("%s  %9.5Lf %9.5Lf %9.5Lf\n", "N", N[3], N[4], N[5]);
   printf("%s  %9.5Lf %9.5Lf %9.5Lf\n", "N", N[6], N[7], N[8]);
#endif

/* Rotation */
   ut1 = SECDAY * (mjd-(int)mjd) + dUT1;	// [sec]
   Tp = (floor(mjd) - 51544.5) / 36525.0;	// [century]
   r = 1.002737909350795 + Tp*(5.9006E-11 - Tp*5.9e-15);
// http://www2.arnes.si/~gljsentvid10/sidereal.htm

//   GMST0 = (6 + 41.0/60 + (50.54841 + Tp*(8640184.812866 +
//                           Tp*(0.093104 - Tp*6.2E-6)))/3600) / 24 * 2*pi;

   GMST0 = 24110.54841 + Tp*(8640184.812866 + Tp*(0.093104 - Tp*6.2e-6)); 
   GMST = GMST0 + r * ut1;			// [sec]
   GMST = fmod(GMST*2*pi/SECDAY, 2*pi);	// [rad]
   GAST = GMST + dpsi*COS(eps) + 
          (0.00264*SIN(Omega) +0.000063*SIN(2*Omega))*secrad;

   rotz(-GAST, R);

#ifdef DEBUG1
   printf("ut1= %8.3Lf  Tp= %10.3f  GMST0= %8.3Lf  GMST= %8.3Lf GAST= %8.3Lf %8.3Lf\n", 
	  ut1, Tp, GMST0*45/ATAN(1.0), GMST*45/ATAN(1.0), GAST*45/ATAN(1.0), 
	  (GAST-GMST)*3600*45/ATAN(1.0));
   printf("%s  %9.5Lf %9.5Lf %9.5Lf\n", "R", R[0], R[1], R[2]);
   printf("%s  %9.5Lf %9.5Lf %9.5Lf\n", "R", R[3], R[4], R[5]);
   printf("%s  %9.5Lf %9.5Lf %9.5Lf\n", "R", R[6], R[7], R[8]);
#endif

   matmult(P, N, PNR);
   matmult(PNR, R, PNR);

#ifdef DEBUG1
   printf("%s  %9.5Lf %9.5Lf %9.5Lf\n", "PNR", PNR[0], PNR[1], PNR[2]);
   printf("%s  %9.5Lf %9.5Lf %9.5Lf\n", "PNR", PNR[3], PNR[4], PNR[5]);
   printf("%s  %9.5Lf %9.5Lf %9.5Lf\n", "PNR", PNR[6], PNR[7], PNR[8]);
#endif

// FIXME: the pole wobble at 0.5" translates to 0.1" at GEO.
//    W = R1(yp)*R2(xp)
//    where xp and yp are the polar coordinates of the celestial ephemeris pole,
//       taken from the IERS bulletin. (NB in the bulletin they are in units
//       of arcseconds, and they must be converted to radians by multiplying
//       by pi/180/3600.)

}
/* Return the Earth's obliquity at mjd, in deg */
/* 1992 Astron Almanac, p. B18, dropping the 2 masec cubic term */
REAL obliquity(REAL mjd)
{
   REAL T;
   T = (mjd - 51544.5) / 36525;  /* centuries since J2000 */
   return(23.439291 + T * (-0.0130042 - 0.00000016 * T));
}


/* Astronomical Almanac (p. D46 in 1992 version) */
/* Low precision moon position from center of Earth: equatorial coords */
void lp_xyzmoon(REAL mjd, VEC *moon)
{
   REAL T, lambda, beta, pie, distance;
   REAL deg_in_radian=45/ATAN(1.0);
   T = (mjd - 51544.5) / 36525.;  /* jul cent. since J2000.0 */

   lambda = 218.32 + 481267.883 * T 
      + 6.29 * SIN((134.9 + 477198.85 * T) / deg_in_radian)
      - 1.27 * SIN((259.2 - 413335.38 * T) / deg_in_radian)
      + 0.66 * SIN((235.7 + 890534.23 * T) / deg_in_radian)
      + 0.21 * SIN((269.9 + 954397.70 * T) / deg_in_radian)
      - 0.19 * SIN((357.5 + 35999.05 * T) / deg_in_radian)
      - 0.11 * SIN((186.6 + 966404.05 * T) / deg_in_radian);
   lambda = lambda / deg_in_radian;
   beta = 5.13 * SIN((93.3 + 483202.03 * T) / deg_in_radian)
      + 0.28 * SIN((228.2 + 960400.87 * T) / deg_in_radian)
      - 0.28 * SIN((318.3 + 6003.18 * T) / deg_in_radian)
      - 0.17 * SIN((217.6 - 407332.20 * T) / deg_in_radian);
   beta = beta / deg_in_radian;
   pie = 0.9508 
      + 0.0518 * COS((134.9 + 477198.85 * T) / deg_in_radian)
      + 0.0095 * COS((259.2 - 413335.38 * T) / deg_in_radian)
      + 0.0078 * COS((235.7 + 890534.23 * T) / deg_in_radian)
      + 0.0028 * COS((269.9 + 954397.70 * T) / deg_in_radian);
   distance = WGS84_a / SIN(pie/deg_in_radian);

// Earth equatorial coords
   moon->x = distance * COS(beta) * COS(lambda);
   moon->y = distance * (0.9175 * COS(beta) * SIN(lambda) - 0.3978 * SIN(beta));
   moon->z = distance * (0.3978 * COS(beta) * SIN(lambda) + 0.9175 * SIN(beta));
}


/* Low precision moon position from center of Earth: ecliptic coordinates */
void lp_xyzmoon_ecliptic(REAL mjd, VEC *moon)
{
   double C, S, incl;

/* Equatorial coordinates */
   lp_xyzmoon(mjd, moon);

/* Rotate to ecliptic coordinates */
   incl = obliquity(mjd) * atan(1)/45;
   C = moon->y;
   S = moon->z;
   moon->y =  COS(incl) * C + SIN(incl) * S;
   moon->z = -SIN(incl) * C + COS(incl) * S;
   return;
}


/* More accurate moon position from center of Earth */
/* from Jean Meeus' *Astronomical Formulae For Calculators */
void xyzmoon(REAL mjd, VEC *moon)
{       
   REAL DEGRAD=ATAN(1.0)/45;
   double pie, dist;  /* horiz parallax */
   double Lpr,M,Mpr,D,F,Om,T,Tsq,Tcb;
   double e,lambda,B,beta,om1,om2;
   double sinx, l, m, n, incl;

/* Approximate correction to from UT to dynamical time, T since 1900 */
   T = (mjd + (33.15+(2.164e-3)*(mjd-36934/*=1960*/))/SECDAY - 15019.5) / 36525.;
   Tsq = T * T;
   Tcb = Tsq * T;

   Lpr = 270.434164 + 481267.8831 * T - 0.001133 * Tsq 
      + 0.0000019 * Tcb;
   M = 358.475833 + 35999.0498*T - 0.000150*Tsq
      - 0.0000033*Tcb;
   Mpr = 296.104608 + 477198.8491*T + 0.009192*Tsq 
      + 0.0000144*Tcb;
   D = 350.737486 + 445267.1142*T - 0.001436 * Tsq
      + 0.0000019*Tcb;
   F = 11.250889 + 483202.0251*T -0.003211 * Tsq 
      - 0.0000003*Tcb;
   Om = 259.183275 - 1934.1420*T + 0.002078*Tsq 
      + 0.0000022*Tcb;

   Lpr = fmod(Lpr, 360.0);
   Mpr = fmod(Mpr, 360.0);
   M = fmod(M, 360.0);
   D = fmod(D, 360.0);
   F = fmod(F, 360.0);
   Om = fmod(Om, 360.0);
	
   sinx =  sin((51.2 + 20.2 * T)*DEGRAD);
   Lpr = Lpr + 0.000233 * sinx;
   M = M - 0.001778 * sinx;
   Mpr = Mpr + 0.000817 * sinx;
   D = D + 0.002011 * sinx;
	
   sinx = 0.003964 * SIN((346.560+132.870*T -0.0091731*Tsq)*DEGRAD);

   Lpr = Lpr + sinx;
   Mpr = Mpr + sinx;
   D = D + sinx;
   F = F + sinx;

   sinx = SIN(Om*DEGRAD);
   Lpr = Lpr + 0.001964 * sinx;
   Mpr = Mpr + 0.002541 * sinx;
   D = D + 0.001964 * sinx;
   F = F - 0.024691 * sinx;
   F = F - 0.004328 * SIN((Om + 275.05 -2.30*T)*DEGRAD);

   e = 1 - 0.002495 * T - 0.00000752 * Tsq;

   M *= DEGRAD;   /* these will all be arguments ... */
   Mpr *= DEGRAD;
   D *= DEGRAD;
   F *= DEGRAD;

   lambda = Lpr + 6.288750 * SIN(Mpr)
      + 1.274018 * SIN(2*D - Mpr)
      + 0.658309 * SIN(2*D)
      + 0.213616 * SIN(2*Mpr)
      - e * 0.185596 * SIN(M) 
      - 0.114336 * SIN(2*F)
      + 0.058793 * SIN(2*D - 2*Mpr)
      + e * 0.057212 * SIN(2*D - M - Mpr)
      + 0.053320 * SIN(2*D + Mpr)
      + e * 0.045874 * SIN(2*D - M)
      + e * 0.041024 * SIN(Mpr - M)
      - 0.034718 * SIN(D)
      - e * 0.030465 * SIN(M+Mpr)
      + 0.015326 * SIN(2*D - 2*F)
      - 0.012528 * SIN(2*F + Mpr)
      - 0.010980 * SIN(2*F - Mpr)
      + 0.010674 * SIN(4*D - Mpr)
      + 0.010034 * SIN(3*Mpr)
      + 0.008548 * SIN(4*D - 2*Mpr)
      - e * 0.007910 * SIN(M - Mpr + 2*D)
      - e * 0.006783 * SIN(2*D + M)
      + 0.005162 * SIN(Mpr - D)
      + e * 0.005000 * SIN(M + D)
      + e * 0.004049 * SIN(Mpr - M + 2*D)
      + 0.003996 * SIN(2*Mpr + 2*D)
      + 0.003862 * SIN(4*D)
      + 0.003665 * SIN(2*D - 3*Mpr)
      + e * 0.002695 * SIN(2*Mpr - M)
      + 0.002602 * SIN(Mpr - 2*F - 2*D)
      + e * 0.002396 * SIN(2*D - M - 2*Mpr)
      - 0.002349 * SIN(Mpr + D)
      + e * e * 0.002249 * SIN(2*D - 2*M)
      - e * 0.002125 * SIN(2*Mpr + M)
      - e * e * 0.002079 * SIN(2*M)
      + e * e * 0.002059 * SIN(2*D - Mpr - 2*M)
      - 0.001773 * SIN(Mpr + 2*D - 2*F)
      - 0.001595 * SIN(2*F + 2*D)
      + e * 0.001220 * SIN(4*D - M - Mpr)
      - 0.001110 * SIN(2*Mpr + 2*F)
      + 0.000892 * SIN(Mpr - 3*D)
      - e * 0.000811 * SIN(M + Mpr + 2*D)
      + e * 0.000761 * SIN(4*D - M - 2*Mpr)
      + e * e * 0.000717 * SIN(Mpr - 2*M)
      + e * e * 0.000704 * SIN(Mpr - 2 * M - 2*D)
      + e * 0.000693 * SIN(M - 2*Mpr + 2*D)
      + e * 0.000598 * SIN(2*D - M - 2*F)
      + 0.000550 * SIN(Mpr + 4*D)
      + 0.000538 * SIN(4*Mpr)
      + e * 0.000521 * SIN(4*D - M)
      + 0.000486 * SIN(2*Mpr - D);
	
   B = 5.128189 * SIN(F)
      + 0.280606 * SIN(Mpr + F)
      + 0.277693 * SIN(Mpr - F)
      + 0.173238 * SIN(2*D - F)
      + 0.055413 * SIN(2*D + F - Mpr)
      + 0.046272 * SIN(2*D - F - Mpr)
      + 0.032573 * SIN(2*D + F)
      + 0.017198 * SIN(2*Mpr + F)
      + 0.009267 * SIN(2*D + Mpr - F)
      + 0.008823 * SIN(2*Mpr - F)
      + e * 0.008247 * SIN(2*D - M - F) 
      + 0.004323 * SIN(2*D - F - 2*Mpr)
      + 0.004200 * SIN(2*D + F + Mpr)
      + e * 0.003372 * SIN(F - M - 2*D)
      + 0.002472 * SIN(2*D + F - M - Mpr)
      + e * 0.002222 * SIN(2*D + F - M)
      + e * 0.002072 * SIN(2*D - F - M - Mpr)
      + e * 0.001877 * SIN(F - M + Mpr)
      + 0.001828 * SIN(4*D - F - Mpr)
      - e * 0.001803 * SIN(F + M)
      - 0.001750 * SIN(3*F)
      + e * 0.001570 * SIN(Mpr - M - F)
      - 0.001487 * SIN(F + D)
      - e * 0.001481 * SIN(F + M + Mpr)
      + e * 0.001417 * SIN(F - M - Mpr)
      + e * 0.001350 * SIN(F - M)
      + 0.001330 * SIN(F - D)
      + 0.001106 * SIN(F + 3*Mpr)
      + 0.001020 * SIN(4*D - F)
      + 0.000833 * SIN(F + 4*D - Mpr)
      + 0.000781 * SIN(Mpr - 3*F)
      + 0.000670 * SIN(F + 4*D - 2*Mpr)
      + 0.000606 * SIN(2*D - 3*F)
      + 0.000597 * SIN(2*D + 2*Mpr - F)
      + e * 0.000492 * SIN(2*D + Mpr - M - F)
      + 0.000450 * SIN(2*Mpr - F - 2*D)
      + 0.000439 * SIN(3*Mpr - F)
      + 0.000423 * SIN(F + 2*D + 2*Mpr)
      + 0.000422 * SIN(2*D - F - 3*Mpr)
      - e * 0.000367 * SIN(M + F + 2*D - Mpr)
      - e * 0.000353 * SIN(M + F + 2*D)
      + 0.000331 * SIN(F + 4*D)
      + e * 0.000317 * SIN(2*D + F - M + Mpr)
      + e * e * 0.000306 * SIN(2*D - 2*M - F)
      - 0.000283 * SIN(Mpr + 3*F);
	
   om1 = 0.0004664 * cos(Om*DEGRAD);        
   om2 = 0.0000754 * cos((Om + 275.05 - 2.30*T)*DEGRAD);
	
   beta = B * (1. - om1 - om2);

   pie = 0.950724 
      + 0.051818 * cos(Mpr)
      + 0.009531 * cos(2*D - Mpr)
      + 0.007843 * cos(2*D)
      + 0.002824 * cos(2*Mpr)
      + 0.000857 * cos(2*D + Mpr)
      + e * 0.000533 * cos(2*D - M)
      + e * 0.000401 * cos(2*D - M - Mpr)
      + e * 0.000320 * cos(Mpr - M)
      - 0.000271 * cos(D)
      - e * 0.000264 * cos(M + Mpr)
      - 0.000198 * cos(2*F - Mpr)
      + 0.000173 * cos(3*Mpr)
      + 0.000167 * cos(4*D - Mpr)
      - e * 0.000111 * cos(M)
      + 0.000103 * cos(4*D - 2*Mpr)
      - 0.000084 * cos(2*Mpr - 2*D)
      - e * 0.000083 * cos(2*D + M)
      + 0.000079 * cos(2*D + 2*Mpr)
      + 0.000072 * cos(4*D)
      + e * 0.000064 * cos(2*D - M + Mpr)
      - e * 0.000063 * cos(2*D + M - Mpr)
      + e * 0.000041 * cos(M + D)
      + e * 0.000035 * cos(2*Mpr - M)
      - 0.000033 * cos(3*Mpr - 2*D)
      - 0.000030 * cos(Mpr + D)
      - 0.000029 * cos(2*F - 2*D)
      - e * 0.000029 * cos(2*Mpr + M)
      + e * e * 0.000026 * cos(2*D - 2*M)
      - 0.000023 * cos(2*F - 2*D + Mpr)
      + e * 0.000019 * cos(4*D - M - Mpr);

   beta *= DEGRAD;
   lambda *= DEGRAD;

/* These are in ecliptic coordinates */
   l = cos(lambda) * cos(beta);    
   m = sin(lambda) * cos(beta);
   n = sin(beta);

   dist = WGS84_a / SIN(pie*DEGRAD);
	
/* Ecliptic coords */
//   moon->x = dist * l;
//   moon->y = dist * m;
//   moon->z = dist * n;

/* Equatorial coords */
   incl = obliquity(mjd) * DEGRAD;
   moon->x = dist * l;
   moon->y = dist * (COS(incl) * m - SIN(incl) * n);
   moon->z = dist * (SIN(incl) * m + COS(incl) * n);
}

/* Low precision formulae for the sun, from Almanac p. C24 (1990) */
void lp_xyzsun(REAL mjd, VEC *sun)
{
   REAL T, L, g, lambda,epsilon;
   REAL deg_in_radian=45/ATAN(1.0);
   T = mjd - 51544.5;
   L = 280.460 + 0.9856474 * T;
   g = (357.528 + 0.9856003 * T)/deg_in_radian;
   lambda = (L + 1.915 * SIN(g) + 0.020 * SIN(2*g)) / deg_in_radian;
   epsilon = (23.439 - 0.0000004 * T) / deg_in_radian;

   sun->x = AU * COS(lambda); 
   sun->y = AU * COS(epsilon) * SIN(lambda); 
   sun->z = AU * SIN(epsilon) * SIN(lambda);
}


/* Convert Kepler params to state vector */
/* From www.cdeagle.com/ccmatlab/ccmatlab.pdf */
void orb2eci(
   double a,	/* Semi-major axis */
   double e,	/* Orbital eccentricity */
   double i,	/* Inclination [rad] */
   double w,	/* (omega) argument of periapsis [rad] */
   double O,	/* (Omega) longitude of ascending node [rad] */
   double n,	/* (nu) *true* anomaly [rad] */
   double gm,	/* GM in desired a,v units */
   double *r,	/* position */
   double *v)	/* velocity */
{
   double p, v0;
/* NB The Matlab doc had this wrong for p... */
   p = a*(1-e*e) / (1+e*cos(n));
   v0 = sqrt(gm/(a*(1-e*e)));

   r[0] = p*(cos(O)*cos(w+n)-sin(O)*cos(i)*sin(w+n));
   r[1] = p*(sin(O)*cos(w+n)+cos(O)*cos(i)*sin(w+n));
   r[2] = p*sin(i)*sin(w+n);
   v[0] = -v0*(cos(O)*(sin(w+n)+e*sin(w)) +
	       sin(O)*cos(i)*(cos(w+n)+e*cos(w)));
   v[1] = -v0*(sin(O)*(sin(w+n)+e*sin(w)) -
	       cos(O)*cos(i)*(cos(w+n)+e*cos(w)));
/* NB The Matlab doc had the sign wrong for vz */
   v[2] =  v0*sin(i)*(cos(w+n)+e*cos(w));
}

/* Mean anomaly to true anomaly */
void ma2ta(
   double e,	/* Eccentricity */
   double M,	/* Mean anomaly [rad]       M = 2*pi*t/P    */
   double *Ecc,	/* Eccentric anomaly [rad]  E = M + e*sin(E)   */
   double *n)	/* (nu) true anomaly [rad]  tan(nu/2) = sqrt((1+e)/(1-e))*tan(E/2)   */
{
   int maxiter=1000, i;
   double f0, f1, f2, f3, D1, D2, dE;
   *Ecc = M + 0.85 * e;
   if(sin(M) < 0) *Ecc = M - 0.85 * e;

   for(i=0; i<maxiter; i++) {
      f0 = *Ecc - e*sin(*Ecc) - M;
      f1 = 1 - e*cos(*Ecc);
      f2 = e*sin(*Ecc);
      f3 = e*cos(*Ecc);
      D1 = -f0 / f1;
      D2 = -f0 / (f1 + D1*f2/2);
      dE = -f0 / (f1 + D1*f2/2 + D2*D2*f3/6);
      *Ecc += dE;
      if(ABS(f0) < 1e-10) break;
   }
   *n = 2*atan2(sqrt(1-e*e)*sin(*Ecc/2), (1-e)*cos(*Ecc/2));
}

/* Return orbital elements of a designated planet k, 1-9 */
/* From ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf */
/* and  ssd.jpl.nasa.gov/?planet_phys_par */
/* NOTE: this approximates Teph as MJD-J2000 so will be off by leap
   seconds since 2000.0 (4sec between 2000 and 2016) */
void planet(int k, double mjd, double *par, double *r, double *v,
	   double *hmag, double *albedo, char *name)
{
   char *planets[9]={"Mercury", "Venus", "Earth", "Mars", "Jupiter",
		     "Saturn", "Uranus", "Neptune", "Pluto"};
// This is the table for 1800-2050; -3000 to +3000 also available
//  a           e            I              L         long.peri.   long.node.
// AU, AU/Cy  rad, rad/Cy  deg, deg/Cy    deg, deg/Cy   deg, deg/Cy  deg, deg/Cy
   double elements[9][6]={
      { 0.38709927,  0.20563593,  7.00497902, 252.25032350,  77.45779628,  48.33076593},
      { 0.72333566,  0.00677672,  3.39467605, 181.97909950, 131.60246718,  76.67984255},
      { 1.00000261,  0.01671123, -0.00001531, 100.46457166, 102.93768193,   0.0},
      { 1.52371034,  0.09339410,  1.84969142,  -4.55343205, -23.94362959,  49.55953891},
      { 5.20288700,  0.04838624,  1.30439695,  34.39644051,  14.72847983, 100.47390909},
      { 9.53667594,  0.05386179,  2.48599187,  49.95424423,  92.59887831, 113.66242448},
      {19.18916464,  0.04725744,  0.77263783, 313.23810451, 170.95427630,  74.01692503},
      {30.06992276,  0.00859048,  1.77004347, -55.12002969,  44.96476227, 131.78422574},
      {39.48211675,  0.24882730, 17.14001206, 238.92903833, 224.06891629, 110.30393684}};
   double derivs[9][6]={
      { 0.00000037,  0.00001906, -0.00594749,149472.67411175,  0.16047689, -0.12534081},
      { 0.00000390, -0.00004107, -0.00078890, 58517.81538729,  0.00268329, -0.27769418},
      { 0.00000562, -0.00004392, -0.01294668, 35999.37244981,  0.32327364,  0.0},
      { 0.00001847,  0.00007882, -0.00813131, 19140.30268499,  0.44441088, -0.29257343},
      {-0.00011607, -0.00013253, -0.00183714,  3034.74612775,  0.21252668,  0.20469106},
      {-0.00125060, -0.00050991,  0.00193609,  1222.49362201, -0.41897216, -0.28867794},
      {-0.00196176, -0.00004397, -0.00242939,   428.48202785,  0.40805281,  0.04240589},
      { 0.00026291,  0.00005105,  0.00035372,   218.45945325, -0.32241464, -0.00508664},
      {-0.00031596,  0.00005170,  0.00004818,   145.20780515, -0.04062942, -0.01183482}};
   double hdata[9]={-0.60, -4.47, -3.86, -1.52, -9.40, -8.88, -7.19, -6.87, -1.0};
   double adata[9]={0.106, 0.65, 0.367, 0.150, 0.52,  0.47, 0.51,  0.41, 0.3};
   double T, a, e, i, O, w, L, periarg, anomean, anoecc, anotrue;
   double dr=atan(1.0)/45, GM=2.95912208e-4/*AU^3/day^2*/;

   T = (mjd - 51544.5) / 36525;

   a = elements[k][0] + T*derivs[k][0];
   e = elements[k][1] + T*derivs[k][1];
   i = elements[k][2] + T*derivs[k][2];
   L = elements[k][3] + T*derivs[k][3];
   w = elements[k][4] + T*derivs[k][4];
   O = elements[k][5] + T*derivs[k][5];

   periarg = w - O;
   anomean = L - w - 360*NINT((L-w)/360);

   par[0] = a;
   par[1] = e;
   par[2] = i;
   par[3] = periarg;
   par[4] = O;
   par[5] = anomean;
   *hmag = hdata[k];
   *albedo = adata[k];
   strcpy(name, planets[k]);

   ma2ta(e, anomean*dr, &anoecc, &anotrue);
   orb2eci(a, e, i*dr, periarg*dr, O*dr, anotrue, GM, r, v);
}


/* Higher? precision formulae for the sun wrt Earth geocenter: equatorial coords */
void xyzsun(REAL mjd, VEC *sun)
{
   double par[6], r[3], v[3], hmag, albedo, incl, DEGRAD=atan(1.0)/45, dmjd=mjd;
   char name[80];
   VEC moon;
   
/* Get the Earth-Moon barycenter's position wrt Sun in ecliptic coords [AU] */
   planet(2, dmjd, par, r, v, &hmag, &albedo, name);

/* Convert to equatorial coords, units of meters, and Earth to Sun */
   incl = obliquity(mjd) * DEGRAD; 
   sun->x = -AU*r[0];
   sun->y = -AU*(COS(incl) * r[1] - SIN(incl) * r[2]);
   sun->z = -AU*(SIN(incl) * r[1] + COS(incl) * r[2]);

/* Get the Moon's position wrt the Earth's geocenter in equatorial coords */
   xyzmoon(mjd, &moon);

/* Sun wrt Earth's geocenter in equatorial coords */
   sun->x += GMMOON/(GMMOON+GMEARTH) * moon.x;
   sun->y += GMMOON/(GMMOON+GMEARTH) * moon.y;
   sun->z += GMMOON/(GMMOON+GMEARTH) * moon.z;
}

/* Higher? precision formulae for the sun wrt Earth geocenter: equatorial coords */
void xyzsun_lpmoon(REAL mjd, VEC *sun)
{
   double par[6], r[3], v[3], hmag, albedo, incl, DEGRAD=atan(1.0)/45, dmjd=mjd;
   char name[80];
   VEC moon;
   
/* Get the Earth-Moon barycenter's position wrt Sun in ecliptic coords [AU] */
   planet(2, dmjd, par, r, v, &hmag, &albedo, name);

/* Convert to equatorial coords, units of meters, and Earth to Sun */
   incl = obliquity(mjd) * DEGRAD; 
   sun->x = -AU*r[0];
   sun->y = -AU*(COS(incl) * r[1] - SIN(incl) * r[2]);
   sun->z = -AU*(SIN(incl) * r[1] + COS(incl) * r[2]);

/* Get the Moon's position wrt the Earth's geocenter in equatorial coords */
   lp_xyzmoon(mjd, &moon);

/* Sun wrt Earth's geocenter in equatorial coords */
   sun->x += GMMOON/(GMMOON+GMEARTH) * moon.x;
   sun->y += GMMOON/(GMMOON+GMEARTH) * moon.y;
   sun->z += GMMOON/(GMMOON+GMEARTH) * moon.z;
}


/* Convert state vector to Keplerian elements and TLE */
// http://space.stackexchange.com/questions/1904/how-to-programmatically-calculate-orbital-elements-using-position-velocity-vecto
// http://celestrak.com/columns/v02n01/
// ADA289281.pdf 

/* Compute Keplerian orbital constants from state vector */
void keparams(ORBIT *orb, REAL GM)
{
   REAL r, x, y, z, L, vr, vt, phix, phiy;
   REAL pi=4*ATAN(1.0);

/* Distance and radial and tangential velocity components */
   r = SQRT(DOT(orb->X,orb->X));
   vr = DOT(orb->V,orb->X) / r;
   x = orb->V.x - vr * orb->X.x/r;
   y = orb->V.y - vr * orb->X.y/r;
   z = orb->V.z - vr * orb->X.z/r;
   vt = SQRT(x*x + y*y + z*z);

/* Angular momentum */
   orb->L.x = XCRS(orb->X,orb->V);
   orb->L.y = YCRS(orb->X,orb->V);
   orb->L.z = ZCRS(orb->X,orb->V);
   L = SQRT(DOT(orb->L,orb->L));

/* Inclination and ascending node */
   orb->incl = ACOS(orb->L.z / L);
   orb->Omega = ATAN2(orb->L.x, -orb->L.y);
   if(orb->Omega < 0) orb->Omega += 2*pi;

/* Pole */
   orb->pra = ATAN2(orb->L.y, orb->L.x);
   if(orb->pra < 0) orb->pra += 2*pi;
   orb->pdec = ASIN(orb->L.z/L);

/* True anomaly = angle CCW from periapse to this point in orbit plane */
   orb->nu = ATAN2(vr, vt-GM/L);
   if(orb->nu < 0) orb->nu += 2*pi;

/* Energy */
   orb->E = 0.5 * DOT(orb->V,orb->V) - GM/r;

/* Semi-major axis and eccentricity */
   orb->a = GM / ABS(2*orb->E);
   orb->ecc = SQRT(1+2*orb->E*L*L/(GM*GM));

/* Period */
   orb->period = 2*pi * SQRT(orb->a*orb->a*orb->a/GM);

/* Periapse from Runge-Lenz vector: |A| = GM*ecc */
   orb->Peri.x = orb->a*(1-orb->ecc)/orb->ecc/GM * (XCRS(orb->V, orb->L)-GM*orb->X.x/r);
   orb->Peri.y = orb->a*(1-orb->ecc)/orb->ecc/GM * (YCRS(orb->V, orb->L)-GM*orb->X.y/r);
   orb->Peri.z = orb->a*(1-orb->ecc)/orb->ecc/GM * (ZCRS(orb->V, orb->L)-GM*orb->X.z/r);

/* Argument of periapse (ascending node to periapse) */
   phix = L * (orb->X.y*orb->L.x - orb->X.x*orb->L.y);
   phiy = (orb->X.x*orb->L.x + orb->X.y*orb->L.y)*orb->L.z - 
           orb->X.z*(orb->L.x*orb->L.x + orb->L.y*orb->L.y);
   orb->omega = ATAN2(phiy, phix) + orb->nu;
   if(orb->omega < 0) orb->omega += 2*pi;
   if(orb->omega > 2*pi) orb->omega -= 2*pi;
}

/* Compute Keplerian elements and derivatives from state vector */
void kepler(ORBIT *orb, REAL GM)
{
   int k;
   REAL r, L, x, y, z, vr, vt;
   REAL pi=4*ATAN(1.0);
   REAL dvr[NDIFEQ], dvt[NDIFEQ], dL[NDIFEQ], phix, phiy;
   REAL dLx[NDIFEQ], dLy[NDIFEQ], dLz[NDIFEQ], dE[NDIFEQ];
   REAL dx[NDIFEQ], dy[NDIFEQ], dz[NDIFEQ];

/* Distance and radial and tangential velocity components */
   r = SQRT(orb->X.x*orb->X.x + orb->X.y*orb->X.y + orb->X.z*orb->X.z);
   vr = (orb->V.x*orb->X.x + orb->V.y*orb->X.y + orb->V.z*orb->X.z) / r;
   x = orb->V.x - vr * orb->X.x/r;
   y = orb->V.y - vr * orb->X.y/r;
   z = orb->V.z - vr * orb->X.z/r;
   vt = SQRT(x*x + y*y + z*z);

/* Angular momentum */
   orb->L.x = orb->X.y*orb->V.z - orb->X.z*orb->V.y;
   orb->L.y = orb->X.z*orb->V.x - orb->X.x*orb->V.z;
   orb->L.z = orb->X.x*orb->V.y - orb->X.y*orb->V.x;
   L = SQRT(orb->L.x*orb->L.x + orb->L.y*orb->L.y + orb->L.z*orb->L.z);

/* Inclination and ascending node */
   orb->incl = ACOS(orb->L.z / L);
   orb->Omega = ATAN2(orb->L.x, -orb->L.y);
   if(orb->Omega < 0) orb->Omega += 2*pi;

/* Pole */
   orb->pra = ATAN2(orb->L.y, orb->L.x);
   if(orb->pra < 0) orb->pra += 2*pi;
   orb->pdec = ASIN(orb->L.z/L);

/* True anomaly = angle CCW from periapse to this point in orbit plane */
   orb->nu = ATAN2(vr, vt-GM/L);
   if(orb->nu < 0) orb->nu += 2*pi;

/* Energy */
   orb->E = 0.5*(orb->V.x*orb->V.x + orb->V.y*orb->V.y + orb->V.z*orb->V.z) -
            GM/r;

/* Semi-major axis and eccentricity */
   orb->a = GM / ABS(2*orb->E);
   orb->ecc = SQRT(1+2*orb->E*L*L/(GM*GM));

/* Period */
   orb->period = 2*pi * SQRT(orb->a*orb->a*orb->a/GM);

/* Periapse (a|1-e|,0,0) in orbit, mix of curposn and forward cross product */
   x = -(orb->X.y*orb->L.z - orb->X.z*orb->L.y) / L;
   y = -(orb->X.z*orb->L.x - orb->X.x*orb->L.z) / L;
   z = -(orb->X.x*orb->L.y - orb->X.y*orb->L.x) / L;
   orb->Peri.x = orb->a/r*ABS(1-orb->ecc) * (orb->X.x*COS(orb->nu)-x*SIN(orb->nu));
   orb->Peri.y = orb->a/r*ABS(1-orb->ecc) * (orb->X.y*COS(orb->nu)-y*SIN(orb->nu));
   orb->Peri.z = orb->a/r*ABS(1-orb->ecc) * (orb->X.z*COS(orb->nu)-z*SIN(orb->nu));

/* Argument of periapse (ascending node to periapse) */
   phix = L * (orb->X.y*orb->L.x - orb->X.x*orb->L.y);
   phiy = (orb->X.x*orb->L.x + orb->X.y*orb->L.y)*orb->L.z - 
           orb->X.z*(orb->L.x*orb->L.x + orb->L.y*orb->L.y);
   orb->omega = ATAN2(phiy, phix) + orb->nu;
   if(orb->omega < 0) orb->omega += 2*pi;
   if(orb->omega > 2*pi) orb->omega -= 2*pi;

/* Derivatives of Keplerian constants of motion wrt state vector */
/* Position */
   dx[0] = 1.0;
   dx[1] = dx[2] = dx[3] = dx[4] = dx[5] = 0.0;
   dy[1] = 1.0;
   dy[0] = dy[2] = dy[3] = dy[4] = dy[5] = 0.0;
   dz[2] = 1.0;
   dz[0] = dz[1] = dz[3] = dz[4] = dz[5] = 0.0;

/* Radial velocity = v.r / r */
   dvr[0] = orb->V.x/r - orb->X.x*vr/(r*r);
   dvr[1] = orb->V.y/r - orb->X.y*vr/(r*r);
   dvr[2] = orb->V.z/r - orb->X.z*vr/(r*r);
   dvr[3] = orb->X.x/r;
   dvr[4] = orb->X.y/r;
   dvr[5] = orb->X.z/r;

/* Tangential speed along orbit = sqrt(v^2-vr^2) */
   dvt[0] = -vr*dvr[0] / vt;
   dvt[1] = -vr*dvr[1] / vt;
   dvt[2] = -vr*dvr[2] / vt;
   dvt[3] = (orb->V.x - vr*dvr[3]) / vt;
   dvt[4] = (orb->V.y - vr*dvr[4]) / vt;
   dvt[5] = (orb->V.z - vr*dvr[5]) / vt;

/* Lx */
   dLx[0] = 0.0;		/* dLx / dx = 0   */
   dLx[1] =  orb->V.z;		/* dLx / dy = vz  */
   dLx[2] = -orb->V.y;		/* dLx / dz = -vy */
   dLx[3] = 0.0;		/* dLx / dvx = 0  */
   dLx[4] = -orb->X.z;		/* dLx / dvy = -z */
   dLx[5] =  orb->X.y;		/* dLx / dvz = y  */
/* Ly */
   dLy[0] = -orb->V.z;		/* dLy / dx = -vz */
   dLy[1] = 0.0;		/* dLy / dy = 0   */
   dLy[2] =  orb->V.x;		/* dLy / dz = vx  */
   dLy[3] =  orb->X.z;		/* dLy / dvx = z  */
   dLy[4] = 0.0;		/* dLy / dvy = 0  */
   dLy[5] = -orb->X.x;		/* dLy / dvz = -x */
/* Lz */
   dLz[0] =  orb->V.y;		/* dLz / dx = vy  */
   dLz[1] = -orb->V.x;		/* dLz / dy = -vx */
   dLz[2] = 0.0;		/* dLz / dz = 0   */
   dLz[3] = -orb->X.y;		/* dLz / dvx = -y */
   dLz[4] =  orb->X.x;		/* dLz / dvy = x  */
   dLz[5] = 0.0;		/* dLz / dvz = 0  */

/* E */
   dE[0] = GM*orb->X.x/(r*r*r);	/* dE / dx = GM x/r^3 */
   dE[1] = GM*orb->X.y/(r*r*r);	/* dE / dy = GM y/r^3 */
   dE[2] = GM*orb->X.z/(r*r*r);	/* dE / dz = GM z/r^3 */
   dE[3] = orb->V.x;		/* dE / dvx = vx */
   dE[4] = orb->V.y;		/* dE / dvy = vy */
   dE[5] = orb->V.z;		/* dE / dvz = vz */

/* Assemble composite derivatives from chain rule */
   for(k=0; k<NDIFEQ; k++) {
      dL[k] = (orb->L.x*dLx[k] + orb->L.y*dLy[k] + orb->L.z*dLz[k]) / L;
/* Inclination */
      orb->deriv[k+0*NDIFEQ] = -(dLz[k] - orb->L.z/L*dL[k]) / 
	 SQRT(L*L - orb->L.z*orb->L.z);
/* Ascending node Omega */
      orb->deriv[k+1*NDIFEQ] = (orb->L.x*dLy[k] - orb->L.y*dLx[k]) /
	 (orb->L.x*orb->L.x + orb->L.y*orb->L.y);
/* Eccentricity */
      orb->deriv[k+2*NDIFEQ] = (L*L*dE[k] + 2*orb->E*L*dL[k]) / 
	 (orb->ecc*GM*GM);
/* Period */
      orb->deriv[k+3*NDIFEQ] = 3*pi * GM * dE[k] /
	 SQRT(-8*orb->E*orb->E*orb->E*orb->E*orb->E);
/* nu = "true anomaly" (not a constant of the motion) */
      orb->deriv[k+5*NDIFEQ] = 
	 ((vt-GM/L)*dvr[k] - vr*(dvt[k]+GM/(L*L)*dL[k])) /
	 (vr*vr+(vt-GM/L)*(vt-GM/L));
/* omega = "argument of periapse" */
      orb->deriv[k+4*NDIFEQ] = orb->deriv[k+5*NDIFEQ] +
	 (phix*(dx[k]*orb->L.x*orb->L.z + orb->X.x*dLx[k]*orb->L.z + orb->X.x*orb->L.x*dLz[k] +
		dy[k]*orb->L.y*orb->L.z + orb->X.y*dLy[k]*orb->L.z + orb->X.y*orb->L.y*dLz[k] - 
		dz[k]*(orb->L.x*orb->L.x + orb->L.y*orb->L.y) -
		2*orb->X.z*(dLx[k]*orb->L.x + dLy[k]*orb->L.y)) -
	  phiy*(dy[k]*orb->L.x*L + orb->X.y*dLx[k]*L + orb->X.y*orb->L.x*dL[k] -
		dx[k]*orb->L.y*L - orb->X.x*dLy[k]*L - orb->X.x*orb->L.y*dL[k]))
	 / (phix*phix+phiy*phiy);
   }

#if 0
   REAL pi=4*ATAN(1.0);
   printf("%8.3LF %8.3LF  %8.3LF %8.3LF  %8.3LF", 
	  orb->pra*180/pi, orb->pdec*180/pi,
	  orb->a*1e-6, orb->ecc, orb->nu*180/pi);
   printf("  %8.3LF %8.3LF %8.3LF\n", 
	  1e-6*orb->Peri.x, 1e-6*orb->Peri.y, 1e-6*orb->Peri.z);
#endif
}


#define TRIPDTACC (1e-10)

/* Calculate orbit and time diffs from three observations and an orbit pole */
int tripdt(REAL pra, REAL pdec, REAL GM, DET *d1, DET *d2, DET *d3,
	   VEC *X,	/* [m] 3D position at d1 */
	   VEC *V,	/* [m/s] 3D velocity at d1 */
	   VEC *P,	/* unit vector of pole */
	   VEC *PERI,	/* unit vector of periapse */
	   REAL *dt2,	/* [s] actual time between d1 and d2 minus orbit time */
	   REAL *dt3)	/* [s] actual time between d1 and d3 minus orbit time */
{
   VEC X1, X2, X3;
   REAL dr=ATAN(1.0)/45, pi=4*ATAN(1.0), rot[9], mat[9];
   REAL r1, r2, r3, dot, l1, l2, l3, dth2, dth3, crs;
   REAL det, a1, a2, a3, a, e, th0, th1;
   REAL b, E1, E2, E3, t1, t2, t3;
   REAL edot, xdot, ydot, x, y;

/* Pole unit vector */
   P->x = COS(pdec) * COS(pra);
   P->y = COS(pdec) * SIN(pra);
   P->z = SIN(pdec);

/* Solve for radii where this plane intercepts each observation */
   dot = d1->X.x * P->x + d1->X.y * P->y + d1->X.z * P->z;
   if(ABS(dot) < TRIPDTACC) {
#ifdef ERR_TRIPDT
      fprintf(stderr, "tripdt: det 1 is perp to pole\n");
#endif
      *dt2 = *dt3 = SECDAY;
      return(-1);
   }
/* l = distance along unit vector from observer to detection */
   l1 = -(d1->O.x * P->x + d1->O.y * P->y + d1->O.z * P->z) / dot;

   dot = d2->X.x * P->x + d2->X.y * P->y + d2->X.z * P->z;
   if(ABS(dot) < TRIPDTACC) {
#ifdef ERR_TRIPDT
      fprintf(stderr, "tripdt: det 2 is perp to pole\n");
#endif
      *dt2 = *dt3 = SECDAY;
      return(-1);
   }
   l2 = -(d2->O.x * P->x + d2->O.y * P->y + d2->O.z * P->z) / dot;

   dot = d3->X.x * P->x + d3->X.y * P->y + d3->X.z * P->z;
   if(ABS(dot) < TRIPDTACC) {
#ifdef ERR_TRIPDT
      fprintf(stderr, "tripdt: det 3 is perp to pole\n");
#endif
      *dt2 = *dt3 = SECDAY;
      return(-1);
   }
   l3 = -(d3->O.x * P->x + d3->O.y * P->y + d3->O.z * P->z) / dot;

   if(l1 < 0 || l2 < 0 || l3 < 0) {
#ifdef ERR_TRIPDT
      fprintf(stderr, "tripdt: plane intersection opposite view direction\n");
#endif
      *dt2 = *dt3 = SECDAY;
      return(-2);
   }

/* 3D locations of the three detections wrt earth center */
   X1.x = d1->O.x + l1 * d1->X.x;
   X1.y = d1->O.y + l1 * d1->X.y;
   X1.z = d1->O.z + l1 * d1->X.z;
/* r = distance from center of the earth */
   r1 = SQRT(X1.x*X1.x + X1.y*X1.y + X1.z*X1.z);

   X2.x = d2->O.x + l2 * d2->X.x;
   X2.y = d2->O.y + l2 * d2->X.y;
   X2.z = d2->O.z + l2 * d2->X.z;
   r2 = SQRT(X2.x*X2.x + X2.y*X2.y + X2.z*X2.z);

   X3.x = d3->O.x + l3 * d3->X.x;
   X3.y = d3->O.y + l3 * d3->X.y;
   X3.z = d3->O.z + l3 * d3->X.z;
   r3 = SQRT(X3.x*X3.x + X3.y*X3.y + X3.z*X3.z);

/* Angles between detection 1 and detections 2 and 3 */
   dth2 = ACOS((X1.x*X2.x + X1.y*X2.y + X1.z*X2.z) / (r1*r2));
   crs = P->x * (X1.y*X2.z - X1.z*X2.y) +
         P->y * (X1.z*X2.x - X1.x*X2.z) +
         P->z * (X1.x*X2.y - X1.y*X2.x);
   if(crs < 0) dth2 *= -1;
   dth3 = ACOS((X1.x*X3.x + X1.y*X3.y + X1.z*X3.z) / (r1*r3));
   crs = P->x * (X1.y*X3.z - X1.z*X3.y) +
         P->y * (X1.z*X3.x - X1.x*X3.z) +
         P->z * (X1.x*X3.y - X1.y*X3.x);
   if(crs < 0) dth3 *= -1;

/* Conic section that includes these three points with focus at earth center */
/* r = a*(1-e^2) / (1+e*cos(th-th0))       th0 = angle from periapse to th0
 * a1 = e*cos(th0)
 * a2 = e*sin(th0)
 * a3 = a*(1-e^2)
 *
 *   ri + ri*a1*cos(thi) + ri*a2*sin(thi) = a3
 *
 *   r1 + r1*a1 = a3
 *   r2*a1*cos(th2) + r2*a2*sin(th2) - r1*a1 = r1 - r2
 *   r3*a1*cos(th3) + r3*a2*sin(th3) - r1*a1 = r1 - r3
 *
 *   (a1) ( r2*cos(th2)-r1    r2*sin(th2) )  = (r1-r2)
 *   (a2) ( r3*cos(th3)-r1    r3*sin(th3) )  = (r1-r3)
 */
   det = ((r2*COS(dth2)-r1) * r3*SIN(dth3) - (r3*COS(dth3)-r1) * r2*SIN(dth2));
   if(ABS(det) < TRIPDTACC*r1) {
#ifdef ERR_TRIPDT
      fprintf(stderr, "tripdt: conic solution is singular\n");
#endif
      *dt2 = *dt3 = SECDAY;
      return(-3);
   }
   a1 = (r3*SIN(dth3) * (r1-r2) - r2*SIN(dth2) * (r1-r3)) / det;
   a2 = (-(r3*COS(dth3)-r1) * (r1-r2) + (r2*COS(dth2)-r1) * (r1-r3)) / det;
   a3 = r1 * (1+a1);

   th0 = ATAN2(a2, a1);
   e = SQRT(a1*a1+a2*a2);
   if(ABS(1-e*e) < TRIPDTACC) {
#ifdef ERR_TRIPDT
      fprintf(stderr, "tripdt: conic solution is a parabola\n");
#endif
      *dt2 = *dt3 = SECDAY;
      return(-4);
   }
   a = a3 / (1-e*e);

/* Evaluate a piece of the conic */
#if 0
   for(i=-10; i<dth2/dr+10; i++) {
      th = i*dr;
      r = a*(1-e*e)/(1+e*COS(th-th0));
      printf("%5d %5.1Lf %10.5Lf %10.5Lf %10.5Lf\n", 
	     i, th/dr, r*1e-6, r*1e-6*COS(th), r*1e-6*SIN(th));
   }
#endif

   if(e > 1) a = -a;	/* Hyperbola */
   if(a < 0) {
#ifdef ERR_TRIPDT
      fprintf(stderr, "tripdt: conic hyperbola has wrong focus\n");
#endif
      *dt2 = *dt3 = SECDAY;
      return(-5);
   }
   b = a * SQRT(ABS(1-e*e));

   if(DTVERB > 0) {
      printf("pole= %6.2Lf %6.2Lf", pra/dr, pdec/dr);
      printf("  l1,2,3= %6.1Lf %6.1Lf %6.1Lf", 1e-6*l1, 1e-6*l2, 1e-6*l3);
      printf("  r1,2,3= %6.1Lf %6.1Lf %6.1Lf", 1e-6*r1, 1e-6*r2, 1e-6*r3);
      printf("  th2,3=  %6.2Lf %6.2Lf", dth2/dr, dth3/dr);
      printf("  a,e,th,b= %7.3Lf %7.3Lf %7.1Lf %7.3Lf\n", 1e-6*a, e, th0/dr, 1e-6*b);
   }
/* Compute the times between 1->2 and 1->3 */
/*
 * Ellipse:
 *   x = a*(cosE - e)
 *   y = b*sinE              b = a*sqrt(1-e^2)
 *   t = sqrt(a^3/GM) * (E - e*sinE)
 * Hyperbola:
 *   x = a*(e-coshE)
 *   y = b*sinhE             b = a*sqrt(e^2-1)
 *   t = sqrt(a^3/GM) * (e*sinhE - E)
 */

   if(e < 1) {		/* Ellipse */
/* eccentric anomalies at the three points */
      E1 = 2*ATAN(SQRT((1-e)/(1+e)) * TAN((0.00-th0)/2));
      E2 = 2*ATAN(SQRT((1-e)/(1+e)) * TAN((dth2-th0)/2));
      E3 = 2*ATAN(SQRT((1-e)/(1+e)) * TAN((dth3-th0)/2));
/* Fix if wrapped across +/-pi (assumes less than 1/2 orbit span) */
      if(ABS(E1-E2)>pi || ABS(E2-E3)>pi || ABS(E3-E1)>pi) {
	 E1 += pi;
	 E2 += pi;
	 E3 += pi;
	 if(E1 > pi) E1 -= 2*pi;
	 if(E2 > pi) E2 -= 2*pi;
	 if(E3 > pi) E3 -= 2*pi;
	 E1 += pi;
	 E2 += pi;
	 E3 += pi;
      }
/* times at the three points relative to periapse */
      t1 = SQRT(a*a*a/GM)*(E1-e*SIN(E1));
      t2 = SQRT(a*a*a/GM)*(E2-e*SIN(E2));
      t3 = SQRT(a*a*a/GM)*(E3-e*SIN(E3));
/* Motion in the orbital plane evaluated at d1 */
      edot = SQRT(GM/(a*a*a)) / (1 - e*COS(E1));
      xdot = -a * SIN(E1) * edot;
      ydot =  b * COS(E1) * edot;

   } else {		/* Hyperbola */
/* eccentric anomalies at the three points */
      E1 = 2*ATANH(SQRT((e-1)/(e+1)) * TAN((0.00-th0)/2));
      E2 = 2*ATANH(SQRT((e-1)/(e+1)) * TAN((dth2-th0)/2));
      E3 = 2*ATANH(SQRT((e-1)/(e+1)) * TAN((dth3-th0)/2));
/* times at the three points relative to periapse */
      t1 = SQRT(a*a*a/GM)*(e*SINH(E1)-E1);
      t2 = SQRT(a*a*a/GM)*(e*SINH(E2)-E2);
      t3 = SQRT(a*a*a/GM)*(e*SINH(E3)-E3);
/* Motion in the orbital plane evaluated at d1 */
      edot = SQRT(GM/(a*a*a)) / (e*COSH(E1) - 1);
      xdot = -a * SINH(E1) * edot;
      ydot =  b * COSH(E1) * edot;
   }
   *dt2 = SECDAY*(d2->obs.mjd - d1->obs.mjd) - (t2 - t1);
   *dt3 = SECDAY*(d3->obs.mjd - d1->obs.mjd) - (t3 - t1);
   if(DTVERB > 0) {
      printf("E1,2,3= %8.3Lf %8.3Lf %8.3Lf", E1,E2,E3);
      printf("   t1,2,3= %8.0Lf %8.0Lf %8.0Lf %8.0Lf %8.0Lf   %8.0Lf %8.0Lf\n", 
	     t1,t2,t3,t2-t1,t3-t1,
	     SECDAY*(d2->obs.mjd-d1->obs.mjd), 
	     SECDAY*(d3->obs.mjd-d1->obs.mjd));
   }

/* Rotation matrix to get vector from world to orbit (periapse at 1,0,0) */
   rotz(pra, mat);
   roty(pi/2-pdec, rot);
   matmult(rot, mat, rot);
   x = rot[0]*X1.x + rot[1]*X1.y + rot[2]*X1.z;
   y = rot[3]*X1.x + rot[4]*X1.y + rot[5]*X1.z;
   th1 = ATAN2(y, x);
   rotz(th1+th0, mat);
   matmult(mat, rot, rot);

/* Periapse (1,0,0) in orbit, apply transpose of rotation */
   PERI->x = rot[0];
   PERI->y = rot[3];
   PERI->z = rot[6];

/* Position at d1 */
   X->x = X1.x;
   X->y = X1.y;
   X->z = X1.z;
/* Velocity at d1: rotate from orbit to world using transpose */
   V->x = rot[0]*xdot + rot[3]*ydot;
   V->y = rot[1]*xdot + rot[4]*ydot;
   V->z = rot[2]*xdot + rot[5]*ydot;

   if(DTVERB > 0) {
      printf("Pole %8.3LF %8.3LF %8.3LF  %8.3LF %8.3LF\n", 
	     P->x,P->y,P->z,pra/dr,pdec/dr);
      printf("Peri %8.3LF %8.3LF %8.3LF\n", PERI->x,PERI->y,PERI->z);
      printf("X    %8.3LF %8.3LF %8.3LF\n", 1e-6*X->x,1e-6*X->y,1e-6*X->z);
      printf("V    %8.3LF %8.3LF %8.3LF\n", 1e-3*V->x,1e-3*V->y,1e-3*V->z);
   }

   return(0);
}
