  AC_PREREQ([2.68])
  AC_INIT([pmtk],[0.1],[pmtk@benjaminfalkner.de])
  AM_INIT_AUTOMAKE([0.1 foreign])
  AC_CONFIG_MACRO_DIR([m4])

  ACLOCAL=" "

  AC_SUBST(NAME,[pmtk])

  AC_CONFIG_SRCDIR([src/config.h.in])
  AC_CONFIG_HEADERS([src/config.h])

  AC_DEFINE([AUTHORS],[""],[Description: AUTHOR])
  AC_DEFINE([PACKAGE],["Protein Motion TK"],[Description: PACKAGE])
  AC_DEFINE([VERSION],["0.1"],[Description: VERSION])
  AC_DEFINE([WWW],[""],[Description: WWW])
  AC_DEFINE([MAIL],["pmtk@benjaminfalkner.de"],[Description: MAIL])
  AC_DEFINE([LICENSE],["GPL3 (GNU General Public License)"],[Description: LICENSE])
  AC_DEFINE([YEARS],["2009-2014"],[Description: YEAR])
  AC_DEFINE([REVISION],["1"],[Description: REVISION])

  AH_BOTTOM([
  	#include <stdio.h>

	#define ERROR_NAME "Error: "
	#define WARNING_NAME "Warning: "
	#define SAY_NAME "  "

	#define FAIL(format, ...) {fprintf(stderr,ERROR_NAME format"(%s:%i)\n", ##__VA_ARGS__,__FILE__,__LINE__);exit(1);} 
	#define WARN(format, ...) {fprintf(stderr,WARNING_NAME format"(%s:%i)\n", ##__VA_ARGS__,__FILE__,__LINE__);}
	#define SAY(format, ...)  {fprintf(stdout,SAY_NAME format"\n", ##__VA_ARGS__);}
	
	#define LICENSE_TEXT \
	"Protein Motion TK (pmtk)  Copyright (C) 2009-2014  Benjamin Falkner \n"\
	"This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.\n"\
	"This is free software, and you are welcome to redistribute it\n"\
	"under certain conditions; type `show c' for details."

	#define HELP_HEAD(PNAME) 	PNAME " - " PACKAGE " (" VERSION "-" REVISION ")\n"	"\n"

	#define USAGE(PNAME,POPTION,PFILE) 	"usage: " PNAME " ["POPTION"] <"PFILE"...>\n"

	#define OPTIONS_HEAD 	"Options:\n"

	#define VERSION_STRING(PNAME) 	PNAME " " PACKAGE " (" VERSION " rev" REVISION ")\n"	"Copyright (C) "YEARS" " AUTHORS "\n"	"Contact: " MAIL"\n"	"Web: " WWW"\n"


	#define str(s) xstr(s)
	#define xstr(s) #s
  ])

  # PLATFORM 
  AC_CANONICAL_HOST
  case "${host_os}" in
  	darwin*) 
		AC_DEFINE([DARWIN],[1],[Description: Make #define Darwin]) dnl
		AC_MSG_NOTICE([Mac OS X])
		;;
	linux*) 
		AC_DEFINE([LINUX],[1],[Description: Make #define Linux]) dnl
		AC_MSG_NOTICE([GNU/LINUX])
		;;
	*BSD*) 
		AC_DEFINE([BSD],[1],[Description: Make #define Bsd]) dnl 
		AC_MSG_NOTICE([BSD (unstable)])
		;;
	*) 
		AC_MSG_ERROR([Your platform is not currently supported]) 
		;;
  esac


  # PROGRAMS
  AC_PROG_CC
  AC_PROG_LN_S
  AC_PROG_AWK
  AC_PROG_CPP
  AC_PROG_INSTALL
  AC_PROG_RANLIB
  AC_PROG_MAKE_SET
  AC_PROG_CXX

  LT_INIT(dlopen shared pic-only)

  #special programs
  AC_CHECK_PROGS(TAR,[tar],none)

  AC_CHECK_PROGS(DOWNLOADER,[wget curl],none)
  case "$ac_cv_prog_DOWNLOADER" in
	wget) AC_SUBST(DOWNLOADER_ARGS,["-O"]);;
	curl) AC_SUBST(DOWNLOADER_ARGS,["-sLO -o"]);;
        *) AC_MSG_WARN([No program for data transfer found]);;
  esac

  AC_PATH_PROG(CHIMERA,chimera,[no],[/Applications/Chimera.app/Contents/MacOS])
  AM_CONDITIONAL([BUILD_CHIMERA], test !"x$ac_cv_path_CHIMERA" = xno)

  # THREADING
  AC_OPENMP

  # HEAEDRS
  AC_CHECK_HEADERS([unistd.h fcntl.h stdlib.h string.h])

  AC_C_INLINE
  AC_TYPE_SIZE_T
  AC_TYPE_SSIZE_T

  AC_FUNC_MALLOC
  AC_FUNC_REALLOC
  AC_CHECK_FUNCS([memset sqrt])


  # MODULES
  AC_SEARCH_LIBS(ssyev_,[mkl atlas lapack],found_lapack=yes,found_lapack=no)
  if test "x$found_lapack" = xno; then
      AC_MSG_ERROR("lapack routines not found")
  fi

  AC_SEARCH_LIBS(read_xtc,[xtc trr xdr],found_xtc=yes,found_xtc=no)
  if test "x$found_xtc" = xno; then
      AC_MSG_WARN([xtc routines not found])
      AC_MSG_NOTICE([make will try to download and build missing library])
  fi
  AM_CONDITIONAL([BUILD_XDR], test "x$found_xtc" = xno)

  AC_SEARCH_LIBS([fftwf_plan_dft_r2c_3d],[fftw3f],found_fftw3=yes,found_fftw3=no)
  AM_CONDITIONAL(HAS_FFTW3, test "x$found_fftw3" = xyes)
  if test "x$found_fftw3" = xno; then
  	AC_MSG_WARN("No FFTW3: - Removing phase space operations\nIf you are using OS X, please check for /usr/local/lib in LD_LIBRARY_PATH")
  else
  	AC_DEFINE([USE_FFTW3],[1],[Description: Fast Fourier Transform])
  fi



  # PYTHON
  AM_PATH_PYTHON([2.4])

  # Adding macro
  AC_DEFUN([AM_CHECK_PYTHON_HEADERS],
 	[AC_REQUIRE([AM_PATH_PYTHON])
 	AC_MSG_CHECKING(for headers required to compile python extensions)
	dnl deduce PYTHON_INCLUDES
	py_prefix=`$PYTHON -c "import sys; print sys.prefix"`
	py_exec_prefix=`$PYTHON -c "import sys; print sys.exec_prefix"`
	PYTHON_INCLUDES="-I${py_prefix}/include/python${PYTHON_VERSION}"
	if test "$py_prefix" != "$py_exec_prefix"; then
	  	PYTHON_INCLUDES="$PYTHON_INCLUDES -I${py_exec_prefix}/include/python${PYTHON_VERSION}"
	fi
	AC_SUBST(PYTHON_INCLUDES)
	dnl check if the headers exist:
	save_CPPFLAGS="$CPPFLAGS"
	CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"
	AC_TRY_CPP([#include <Python.h>],dnl
	   	[AC_MSG_RESULT(found)
	   	$1],dnl
	   	[AC_MSG_RESULT(not found)
	   	$2])
	   	CPPFLAGS="$save_CPPFLAGS"
	   	])
  AM_CHECK_PYTHON_HEADERS(,[AC_MSG_ERROR(could not find Python headers)])

  # Debug
  AC_ARG_ENABLE(debug,AC_HELP_STRING(--enable-debug, create a debug build),found_debug=$enable_debug)
  AM_CONDITIONAL(IS_DEBUG, test "x$found_debug" = xyes)

  # Static
  # dnlAC_ARG_ENABLE(static,AC_HELP_STRING(--enable-static, create a static build),[LDFLAGS="$LDFLAGS -static"])



  AC_CONFIG_FILES([Makefile src/Makefile src/lib/Makefile python/Makefile xdr/Makefile chimera/Makefile]) #Makefiles
  AC_OUTPUT
