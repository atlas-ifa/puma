CC= cc
# CFLAGS= -g -Wall
CFLAGS= -O3 -g -Wall

SORTDIR=../sortlib
KDDIR=../kdtree

# Use difsys.c for libpuma?
# DIFSYS=difsys.o
DIFSYS=
ifneq ($(DIFSYS),)
CFLAGS+= -DUSE_DIFSYS
endif


all: puma mpclat pumalink

libpuma.a: puma.h orbit.h libpuma.o orbutil.o pumamrq.o $(DIFSYS)
	ar rv libpuma.a libpuma.o orbutil.o pumamrq.o $(DIFSYS)

puma: puma.o libpuma.a
	$(CC) $(CFLAGS) -L$(SORTDIR) puma.o libpuma.a -ltsort -lm -o puma

pumalink.o: pumalink.c pumalink.h
	$(CC) $(CFLAGS) -I$(KDDIR) -c pumalink.c

pumalink: pumalink.o libpuma.a
	$(CC) $(CFLAGS) -L$(SORTDIR) -L$(KDDIR) pumalink.o libpuma.a $(KDDIR)/libkdtree.a -ltsort -lm -o pumalink

mpclat: mpclat.o
	$(CC) $(CFLAGS) mpclat.o -lm -o mpclat

libpuma.o: orbit.h puma.h
puma.o: orbit.h

clean:
	rm *.o puma libpuma.a mpclat pumalink

install:
	install -p -m 0775 puma $(ATLAS_HOME)/bin
	install -p -m 0775 mpclat $(ATLAS_HOME)/bin
	install -m 0644 puma.man $(ATLAS_HOME)/man/man1/puma.1
	install -p -m 0644 puma.h $(ATLAS_HOME)/include
