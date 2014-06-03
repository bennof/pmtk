SUBDIRS = lib

# Global Variables
ACLOCAL_AMFLAGS  =
AM_CFLAGS        = --pedantic -Wall -std=gnu99 -O2 -Ilib
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
AM_LDFLAGS += $(shell make -C ../xdr getstaticlibpath) 
endif

AM_LDFLAGS += lib/libpmtk.a
bin_PROGRAMS = pmdfilter pmdlintraj pmdpca pmmconvert pmmlintraj pmmpca pmmrmsd

pmdfilter_SOURCES  = pmdfilter.c   
pmdlintraj_SOURCES = pmdlintraj.c
pmdpca_SOURCES     = pmdpca.c
pmmconvert_SOURCES = pmmconvert.c
pmmlintraj_SOURCES = pmmlintraj.c
pmmpca_SOURCES     = pmmpca.c
pmmrmsd_SOURCES    = pmmrmsd.c     

