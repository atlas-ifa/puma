CC= gcc
FC= gfortran
CPP= g++
CFLAGS= -O3 -Wall
CLIBS= -lstdc++

libtsort.a: tsort.o
	ar r libtsort.a tsort.o

tsort.o: tsort.cpp tsort.h
	$(CPP) $(CFLAGS) -c -o tsort.o tsort.cpp

install:
	install -p -m 0775 libtsort.a $(ATLAS_HOME)/lib/
	install -p -m 0775 tsort.h $(ATLAS_HOME)/include/
	install -p -m 0664 tsort.man $(ATLAS_HOME)/man/man1/tsort.1

all: tsort

clean:
	rm -f *.o *.a


