#!/bin/bash
# Parse MPC80 format and write TRD9
# Syntax: mpc80_to_trd9.sh [options] < MPC > TRD9

# INSTRUCTIONS:
# -------------
# 
# The TRD9 format required by pumalink has 9 fields (anything beyond is ignored)
# 
#     1. MJD
#     2. J2000 RA [deg]
#     3. J2000 Dec [deg]
#     4. cross-track astrometric uncertainty [arcsec]
#     5. along-track astrometric uncertainty [arcsec]
#     6. observatory east longitude [deg]
#     7. observatory latitude [deg]
#     8. observatory elevation above the WGS80 geoid [m]
#     9. unique ID, no white space, no commas
# 
# This script translates a file traditional MPC80 format lines, one per
# detection, into TRD9 format.  For the astrometric errors it can use
# a value that the user provides (0.2" by default) for all detections or
# it can use an "astrometric error model" described below.
#
# For the ID it assembles a string called "epE.NNNN" (for example ep1.0007),
# where E is an integer that is provided by the user (default 1) and NNNN
# is a running count of the MPC80 lines, starting at 0001.
# 
# Everything after the 9 fields will be ignored by the TRD9 reader, but
# the script provides a ampersand as a visual separator, the magnitude from
# the MPC80 file, the object's name provided in the MPC80, the bandpass for
# the magnitude, the reporting MPC site, and an asterisk for discoveries.
# 
# A constant astrometric error for all detections is probably not
# satisfactory because in fact pumalink really does depend on decent
# astrometric errors.  Options for the user:
# 
# 1. Accept a constant astrometric error, possibly provided as 
# 
#   mpc80_to_trd90.sh xerr=0.1 < file.mpc80 > file.trd9
# 
# 2. Accept a constant astrometric error, but then process the output through
#    some other filter that can amend the errors on a per-detection basis
# 
# 3. Exploit the mpc80_to_trd9.sh "astrometric model" that requires input
#    of a limiting magnitude, PSF FWHM [arcsec], and error floor, e.g.
# 
#   mpc80_to_trd90.sh xerr=0.1 mim=21 fwhm=2.5 < file.mpc80 > file.trd9
# 
#    This will cause mpc80_to_trd90.sh to calculate an astrometric error
#    on a per-detection basis as
#
#   error["] = fwhm["] * sqrt(0.04*exp(0.4*log(10)*(m-mlim)) + xerr*xerr) }
#
#    where 0.04*exp(0.4*log(10)*(m-mlim)) is just the square of the
#    magnitude uncertainty, assuming it is 0.2 at mlim.
#
#    Obviously this model will not apply to a heterogeneous set of MPC80
#    sites, so you need to use it individually, and if you do that you
#    need to be careful to keep the TRD9 IDs unique.
#
# A final choice for the user is to provide an observatory site code,
# east longitude [deg], latitude [deg], and elevation [m].  This avoids
# having to access the MPC sitecode URL, and the sitecode provided should
# match the code in the MPC80 line.
#
# For any input MPC80 line that does not match an MPC site code
# mpc80_to_trd9.sh will supply 0 for longitude, latitude, and elevation.
# 
# SUMMARY OF USEFUL OPTIONS:
# --------------------------
# 
# 1. Provide an epoch number to distinguish one group of IDs from another
# 
#     mpc80_to_trd90.sh epoch=3 < file.mpc80 > file.trd9
# 
# 2. Provide a fixed astrometric error [arcsec]
# 
#     mpc80_to_trd90.sh xerr=0.1 < file.mpc80 > file.trd9
# 
# 3. Provide an astrometric error model
# 
#     mpc80_to_trd90.sh xerr=0.1 mlim=21 fwhm=2.5 < file.mpc80 > file.trd9
# 
# 4. Provide an observatory specification for LAX (to match MPC80 input data)
# 
#     mpc80_to_trd90.sh sitecode=LAX lng=-118.4087 lat=33.9467 elev=20 < file.mpc80 > file.trd9
# 

# Site defaults that can be set instead of MPC lookup (not yet implemented)
sitecode=""
lng=0
lat=0
elev=0

# Astrometric error [arcsec]
xerr=0.2

# PSF full width half maximum [arcsec] for astrometric error model
fwhm=0
# 5-sigma limiting magnitude for astrometric error model
mlim=0

# What epoch label should we use?
epoch=1

# Where can observatory sitecodes be found?
MPC_site_url="https://minorplanetcenter.net/iau/lists/ObsCodes.html"

VERBOSE=0
CLEAN=2

# Get arguments
eval $@

if [[ $VERBOSE -gt 1 ]] ; then set -x ; fi

# Temporary file to hold sitecode,lng,lat,elev
mpclle="$(mktemp /tmp/mpc_site.XXXXXXXXXX)"

# Skip the MPC lookup and use sitecode lng lat elev from the command line?
if [[ -n $sitecode ]] ; then
  echo $sitecode $lng $lat $elev > $mpclle

# Download the current list of all MPC sites from $MPC_site_url
else
  mpcode="$(mktemp /tmp/mpc_site.XXXXXXXXXX)"

# Get the obs code list
  wget -O- -q --no-check-certificate $MPC_site_url > "$mpcode"

# Convert them all to lng, lat, elev
  cat "$mpcode" | while read mpcsite ; do
    if [[ ${mpcsite:0:1} == "<" ]] ; then continue ; fi
    if [[ ${mpcsite:0:4} == "Code" ]] ; then continue ; fi
    lng=${mpcsite:4:9}
    cos=${mpcsite:13:8}
    sin=${mpcsite:21:9}
    name=${mpcsite:30}
    lathgt="$(mpclat $cos $sin)"
    echo ${mpcsite:0:3} $lng $lathgt $name
  done > $mpclle
fi

# Process all the MPC lines into 
#  MJD  RA  Dec  err  err  obslng  obslat  elev  ID &  m    Object  f site
# if its MPC code exists
awk -v tz=$(date +%z) -v lle=$mpclle -v err=$xerr -v mlim=$mlim -v fw=$fwhm -v ep=$epoch '
# BEGIN by reading the entire mpclle file into arrays for site code lookup
    BEGIN{if(length(lle)>0) { while(( getline line<lle) > 0 ) {
               split(line, a, " ");
               lng[a[1]]=a[2]<180?a[2]:a[2]-360;
               lat[a[1]]=a[3]; elev[a[1]]=a[4]}   }
          idnum=0;     # initialize ID number
         }
     {idnum++;         # increment ID number

# Assemble the MJD
      date=substr($0,16,17);   # pluck out the date and time field
      Y=substr(date,1,4); M=substr(date,6,2); D=substr(date,9,2);
      h=int(24*substr(date,11)); m=int(60*(24*substr(date,11)-h)); 
      s=60*(60*(24*substr(date,11)-h)-m); fracsec=s-int(s);
      tzsec=substr(tz,2,2)*3600+substr(tz,4,2)*60;
      if(substr(tz,1,1)=="-") tzsec*=-1;
      datespec=sprintf("%04d %02d %02d %02d %02d %02d", Y,M,D,h,m,int(s));
      mjd=(mktime(datespec)+fracsec+tzsec)/86400.0+40587.0;

# Assemble the RA and Dec from the various fields
      ra=substr($0,33,12); dsgn=substr($0,45,1); dec=substr($0,46,12); 
      r=15*(substr(ra,1,2)+substr(ra,4,2)/60+substr(ra,7)/3600);
      d=(substr(dec,1,2)+substr(dec,4,2)/60+substr(dec,7)/3600);
      if(dsgn == "-") d=-d;

# magnitude field, filter name ("x" for blank filter name)
      m=substr($0,66,5); filt=substr($0,71,1); if(filt==" ") filt="x";

# if both mlim>0 and fw>0 apply an astrometric model instead of constant error
      xerr=terr=err;   # default, just in case
      if(mlim>0 && fw>0 && m>0) {
        dm=0.2*exp(0.2*log(10)*(m-mlim)); xerr=terr=fw*sqrt(dm*dm+err*err) }

# Object name, discovery *, site code
      name=substr($0,6,7); star=substr($0,13,1); obscode=substr($0,78,3);
      printf "%12.6f %9.5f %9.5f %6.3f %6.3f %10.5f %9.5f %5.0f ep%d.%04d & %5.2f %s %s %s%s\n", 
             mjd, r, d, xerr, terr, lng[obscode], lat[obscode], elev[obscode], 
             ep, idnum, m, name, filt, obscode, star}'

# Clean up
if [[ $CLEAN -gt 0 ]] ; then
  if [[ mpclookup -gt 0 ]] ; then rm -f $mpcode $mpclle ; fi
fi

exit 0

