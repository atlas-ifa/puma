# PUMA
PUMA, short for "Position Using Motion under Acceleration", is a collection of C and C++ source code to perform fast, accurate short-art orbit determination and linking using algorithms described in Tonry (2023) in prep.

The code consists of two executables:

* puma - fit short arc orbits
* pumalink - evaluate inter-night and possibly inter-observatory linkages of pairs of detections

PUMA's input is a text file of detections in "TRD9" (time, Right Ascension, Declination, cross-track error, along-track error, observatory east longitude, observatory east latitude, observatory elevation, id) format.  We provide some utility Bash scripts to convert MPC80 to TRD format.

There are four directories present:

* puma - the puma, pumalink, and MPC conversion functions
* sortlib - the ATLAS library of sorting functions
* kdtree - a kdtree library from John Tsiombikas
* examples - a set of examples that illustrate the use of puma

These should be compiled in dependency order:

    cd kdtree ; make ; cd ..
    cd sortlib ; make ; cd ..
    cd puma ; make ; cd ..

Man pages for the programs with extension ".man".  It is the user's
responsibility to decide where the executables should live and to
copy them there.

We provide examples for puma and pumalink in the examples directory,
with README files that describe the operation.
