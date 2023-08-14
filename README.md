# PUMA
PUMA, short for "Position Using Motion with Acceleration", is a collection of C and C++ source code to perform fast, accurate short-arc orbit determination and linking using algorithms described in Tonry (2023) in prep.

The code consists of two executables:

* puma - fit short arc orbits
* pumalink - evaluate linkages of pairs of detections, potentially inter-night and inter-observatory 

PUMA's input is a text file of detections in "TRD9" (MJD, right ascension, declination, cross-track error, along-track error, observatory east longitude, observatory latitude, observatory elevation, ID) format.  We provide some utility Bash scripts to convert MPC80 to TRD9 format.

There are four directories present:

* puma - the puma, pumalink, and MPC conversion functions
* sortlib - the ATLAS library of sorting functions
* kdtree - a kdtree library from John Tsiombikas (repo at https://github.com/jtsiomb/kdtree) 
* examples - a set of examples that illustrate the use of puma

These should be compiled in dependency order:

    cd kdtree ; make ; cd ..
    cd sortlib ; make ; cd ..
    cd puma ; make ; cd ..

Man pages for the programs exist with extension ".man".  It is the user's responsibility to decide where the executables should live and to copy them there.

We provide examples for puma and pumalink in the examples directory, with README files that describe the operation.
