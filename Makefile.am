AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4

if BUILD_XDR
OPT_XTC = xdr
endif


if BUILD_CHIMERA
OPT_CHIMERA = chimera
endif

SUBDIRS = $(OPT_XTC) src python $(OPT_CHIMERA)
