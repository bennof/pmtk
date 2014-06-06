SUBDIRS = lib

# Global Variables
ACLOCAL_AMFLAGS  =
AM_CFLAGS        = --pedantic -Wall -std=c99 -O2 -Ilib
AM_LDFLAGS       =  #:w-Wl,--whole-archive lib/libpmtk.a -Wl,--no-whole-archive

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
AM_LDFLAGS += $(shell make -C ../xdr getstaticlibpath) 
endif



#Build binaries
bin_PROGRAMS = pmdfilter pmdlintraj pmdpca pmmconvert pmmlintraj pmmpca pmmrmsd

pmdfilter_SOURCES  = pmdfilter.c   
pmdfilter_LDFLAGS  = lib/libpmtk.a $(AM_LDFLAGS)
pmdlintraj_SOURCES = pmdlintraj.c
pmdlintraj_LDFLAGS  = lib/libpmtk.a $(AM_LDFLAGS)
pmdpca_SOURCES     = pmdpca.c
pmdpca_LDFLAGS  = lib/libpmtk.a $(AM_LDFLAGS)
pmmconvert_SOURCES = pmmconvert.c
pmmconvert_LDFLAGS  = lib/libpmtk.a $(AM_LDFLAGS)
pmmlintraj_SOURCES = pmmlintraj.c
pmmlintraj_LDFLAGS  = lib/libpmtk.a $(AM_LDFLAGS)
pmmpca_SOURCES     = pmmpca.c
pmmpca_LDFLAGS  = lib/libpmtk.a $(AM_LDFLAGS)
pmmrmsd_SOURCES    = pmmrmsd.c     
pmmrmsd_LDFLAGS  = lib/libpmtk.a $(AM_LDFLAGS)


#Build python Module
pyexec_LTLIBRARIES = py@NAME@.la
py@NAME@_la_CFLAGS   = $(PYTHON_INCLUDES) $(AM_CFLAGS)
py@NAME@_la_SOURCES  = pypmtk.c
py@NAME@_la_LDFLAGS  = -avoid-version -module  $(AM_LDFLAGS)
py@NAME@_la_LIBADD   = lib/libpmtk.a

