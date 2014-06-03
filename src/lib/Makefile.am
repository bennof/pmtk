# Global Variables
ACLOCAL_AMFLAGS  =
AM_CFLAGS        = --pedantic -Wall -std=gnu99 -O2 -I..
AM_LDFLAGS       = 

# Debug section
if IS_DEBUG
AM_CFLAGS += -g -ggdb -O0
AM_CFLAGS += -Wno-long-long -Wall -W -Wnested-externs -Wformat=2
AM_CFLAGS += -Wmissing-prototypes -Wstrict-prototypes -Wmissing-declarations
AM_CFLAGS += -Wwrite-strings -Wshadow -Wpointer-arith -Wsign-compare
AM_CFLAGS += -Wundef -Wbad-function-cast -Winline -Wcast-align
endif

# Add FFTW to linker
if HAS_FFTW3
AM_LDFLAGS += -lfftw3f
endif

if BUILD_XDR
AM_LDFLAGS += $(shell make -C ../../xdr getstaticlibpath)
endif

noinst_LIBRARIES = lib@NAME@.a
lib@NAME@_a_SOURCES = calc.h	config.h model.c                         mrc.h pdb.h protein.c			 calc.c	calctools.c density.c			 mrc.c pdb.c pmtk.h vec.h
lib@NAME@_a_CFLAGS = $(AM_CFLAGS)

lib_LTLIBRARIES = lib@NAME@.la
lib@NAME@_la_SOURCES = calc.h config.h model.c                          mrc.h pdb.h protein.c			  calc.c calctools.c density.c			  mrc.c pdb.c pmtk.h vec.h
include_HEADERS = pmtk.h vec.h



