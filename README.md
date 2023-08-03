# PUMA
PUMA, short for "Position using Motion and Acceleration", is a collection of C and C++ source code to perform fast, accurate short-arc orbit determination and linking using
algorithms described in Tonry (2023) in prep.

The code consists of two executables:

* puma - fit short arc orbits
* pumalink - evaluate inter-night and possibly inter-observatory linkages of pairs of detections

PUMA's input is a text file of detections in "TRD" (time, Right Ascension, Declination) format.  TRD also includes observer longitude, latitude and elevation.  We provide some utility Bash scripts to convert MPC80 to TRD format.
