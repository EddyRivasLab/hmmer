# Process this file with autoconf to produce the HMMER configure script.
#
# Note that this is synchronized with Easel's configure script,
# so that HMMER configures Easel in addition to itself; this saves
# having to do a separate ./configure in Easel. That is,
# everything that appears in Easel's configure script must also
# appear here. 
#
# reminders to save re-reading autoconf manual for the n'th time:
#   - output variables:
#        are defined here as normal shell variables, e.g. FOO="my string"
#        are made into output variables by calling AC_SUBST(FOO)
#        any occurrence of @FOO@ in an output file is then substituted
#        I think this only happens in files we assign w/ AC_CONFIG_FILES;
#        that's the Makefile.in's.
#
#   - C preprocessor symbols:
#        are defined here by calling AC_DEFINE(FOO) or AC_DEFINE(FOO, [42])
#        then #undef FOO lines in easel.h.in become #define FOO or #define FOO 42
#        I think this only happens in header files that we assign
#        w/ AC_CONFIG_HEADERS -- which means, easel.h.in
#
# SRE, Mon Oct  5 14:55:45 1998
# SVN $Id$
# xref autoconf macro archive: //www.gnu.org/software/ac-archive/
#
# GNU recommends the following order:
#   1. autoconf requirements
#   2. AC_INIT
#   3. info on the package
#   4. checks for programs
#   5. checks for libraries
#   6. checks for header files
#   7. checks for types
#   8. checks for structures
#   9. checks for compiler characteristics
#  10. checks for library functions
#  11. checks for system services
#  12. AC_CONFIG_FILES
#  13. AC_OUTPUT


AC_PREREQ(2.60)
AC_INIT(HMMER, 3.0, eddys@janelia.hhmi.org, hmmer)
AC_MSG_NOTICE([Configuring HMMER for your system.])


################################################################
# 3. Info on the package
#
# The four AC_INIT args set these output variables and preprocessor symbols:
#     PACKAGE_NAME      <package>     e.g. "HMMER"
#     PACKAGE_VERSION   <version>     e.g. "3.0"
#     PACKAGE_BUGREPORT <bug-report>  e.g. "eddys@janelia.hhmi.org"
#     PACKAGE_TARNAME   <tarname>     e.g. "hmmer"
# From them, it derives one more:
#     PACKAGE_STRING    <package> <version>, e.g. "HMMER 3.0"
# We also define additional variables:
#     HMMER_RELCODE     unique CVS tag without .'s: e.g. "hmmer3_0"
#     HMMER_DATE        release date: e.g. "January 2007"
#     HMMER_COPYRIGHT   one-line copyright string
#     HMMER_LICENSE     one-line license string
#     HMMER_LICENSETAG  which license to bundle from Licenses/ subdirectory.
#     HMMER_VERSION     copy of version code, e.g. "3.0"
#     HMMER_URL         URL home for HMMER.
# And we have to define the relevant package variables for Easel as well.
#
# We avoid using AC_INIT's PACKAGE_ variables anywhere, because we want to be able
# to use HMMER as a library inside other packages, with no name clashes.
################################################################

HMMER_RELCODE="hmmer3_0"
HMMER_DATE="June 2007"
HMMER_COPYRIGHT="Copyright (C) 1992-2007 HHMI Janelia Farm Research Campus"
HMMER_LICENSE="Freely distributed under the GNU General Public License (GPL)"
HMMER_LICENSETAG=gnu
HMMER_VERSION=$PACKAGE_VERSION
HMMER_URL="http://hmmer.janelia.org/"

EASEL_RELCODE="easel0_1"
EASEL_DATE="June 2007"
EASEL_COPYRIGHT="Copyright (C) 2004-2007 HHMI Janelia Farm Research Campus"
EASEL_LICENSE="Freely distributed under the Janelia Software License."
EASEL_LICENSETAG=jsl
EASEL_VERSION="0.1"
EASEL_URL="http://easel.janelia.org/"

# Output variables  (AC_OUTPUT replaces @var@ in input files, such as Makefiles)
AC_SUBST(HMMER_RELCODE)
AC_SUBST(HMMER_DATE)
AC_SUBST(HMMER_COPYRIGHT)
AC_SUBST(HMMER_LICENSE)
AC_SUBST(HMMER_LICENSETAG)
AC_SUBST(HMMER_VERSION)
AC_SUBST(HMMER_URL)

AC_SUBST(EASEL_RELCODE)
AC_SUBST(EASEL_DATE)
AC_SUBST(EASEL_COPYRIGHT)
AC_SUBST(EASEL_LICENSE)
AC_SUBST(EASEL_LICENSETAG)
AC_SUBST(EASEL_VERSION)
AC_SUBST(EASEL_URL)

# Preprocessor symbols (replace #undefs in p7config.h)
AC_DEFINE_UNQUOTED(HMMER_DATE,      "$HMMER_DATE")
AC_DEFINE_UNQUOTED(HMMER_COPYRIGHT, "$HMMER_COPYRIGHT")
AC_DEFINE_UNQUOTED(HMMER_LICENSE,   "$HMMER_LICENSE")
AC_DEFINE_UNQUOTED(HMMER_VERSION,   "$HMMER_VERSION")

AC_DEFINE_UNQUOTED(EASEL_DATE,      "$EASEL_DATE")
AC_DEFINE_UNQUOTED(EASEL_COPYRIGHT, "$EASEL_COPYRIGHT")
AC_DEFINE_UNQUOTED(EASEL_LICENSE,   "$EASEL_LICENSE")
AC_DEFINE_UNQUOTED(EASEL_VERSION,   "$EASEL_VERSION")

AC_DEFINE(eslLIBRARY)

# Figure out what host we're compiling on.
# Three GNU scripts must be included in the distro: 
#       install.sh, config.guess, config.sub
# This sets four shell variables:
#       host            example: i686-pc-linux-gnu      
#       host_cpu        example: i686
#       host_vendor     example: pc
#       host_os         example: linux-gnu
#
# These are used later in the configure by the ACX_PTHREAD macro to determine
# system-specific details of threads libraries.
#
AC_CANONICAL_HOST

# Check if we're a binary distribution. 
# Trigger is existence of binaries/, nonexistence of src/.
# 
# For binary distro, we do minimal configuration, setting version information 
# and installation path names, making ./configure independent of having a working
# C compiler.
#
# Implements a sneaky if statement that wraps most of the configure.ac M4 code.
#
if ((! test -d binaries) && test -d src); then 
AC_MSG_NOTICE([Full HMMER source distribution - using full configuration])

# Choose the compiler
#
AC_PROG_CC
if test "$CC" = xlc ; then	# IBM's xlc?
  CFLAGS=$CFLAGS" -O3 -qarch=auto -qtune=auto"
fi

################################################################
# Process the ./configure command line
################################################################

# --enable-debugging=x    - set debugging level to <x> (1-3)
#
# At all levels, including 0, replaces CFLAGS w/ "-g -Wall" (so it assumes gcc).
# Sets the p7_DEBUGLEVEL preprocessor symbol to <x>
#
AC_ARG_ENABLE(debugging,
[  --enable-debugging      set CFLAGS for debugging
  --enable-debugging=x    also set diagnostics level to <x> (1-3) (3 = most verbose)],
[ case $enable_debugging in
   yes)  AC_MSG_NOTICE([enabled debugging diagnostics level 0 (CFLAGS only, no verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(p7_DEBUGLEVEL, 0)
	 AC_DEFINE(eslDEBUGLEVEL, 0)
         ;;
     1)  AC_MSG_NOTICE([enabled debugging diagnostics level 1 (low verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(p7_DEBUGLEVEL, 1)
	 AC_DEFINE(eslDEBUGLEVEL, 1)
         ;;
     2)  AC_MSG_NOTICE([enabled debugging diagnostics level 2 (moderate verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(p7_DEBUGLEVEL, 2)
	 AC_DEFINE(eslDEBUGLEVEL, 2)
         ;;
     3)  AC_MSG_NOTICE([enabled debugging diagnostics level 3 (high verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(p7_DEBUGLEVEL, 3)
	 AC_DEFINE(eslDEBUGLEVEL, 3)
         ;;
    no)  AC_MSG_NOTICE([debugging diagnostics disabled])
         AC_DEFINE(p7_DEBUGLEVEL, 0)
	 AC_DEFINE(eslDEBUGLEVEL, 0)
         ;;
     *)  echo "Ignoring unknown argument to --enable-debugging: $enable_debugging"
         ;;
esac])

# --enable-gcov      - compile for code coverage testing
#
# Replaces CC with "gcc" and CFLAGS with "-g -Wall -fprofile-arcs -ftest-coverage".
# Running programs (notably the test suite) will then generate .gcda files,
# which the gcov tool can read.
#
AC_ARG_ENABLE(gcov, 
  [  --enable-gcov           Compile with instrumentation for code coverage testing], 
  [ 
    AC_MSG_NOTICE([Compiling with gcov instrumentation.])
    CC="gcc"
    CFLAGS="-g -Wall -fprofile-arcs -ftest-coverage"
  ])


# --enable-impl      - choose an optimized implementation
#
impl_choice="fast"
impl_cflags=""
AC_ARG_ENABLE(impl,
              [AS_HELP_STRING([--enable-impl],
                [enable an optimized DP implementation (default=fast)]) ],
              [ $impl_choice=$enable_impl])

case $impl_choice in
reference) AC_MSG_NOTICE([Using reference (unoptimized) DP implementation])
           AC_DEFINE(p7_IMPL_REFERENCE)
           ;;

fast)      AC_MSG_NOTICE([Using standard (optimized) DP implementation])
           AC_DEFINE(p7_IMPL_FAST)
           ;;

altivec)   AC_MSG_NOTICE([Using Altivec DP implementation for PowerPC (Erik Lindahl)])
           AC_DEFINE(p7_IMPL_ALTIVEC)
	   if test "$CC" != xlc; then 
	     if test "$ac_cv_build" = powerpc64-unknown-linux-gnu; then 
  	       impl_cflags="-maltivec -mabi=altivec"      # IBM gcc
  	     else
	       impl_cflags="-faltivec -mcpu=G5 -mtune=G5" # Apple cc
	     fi
 	   else
	     impl_cflags="-qaltivec"  # IBM xlc
 	   fi
	   ;;

buhler)    AC_MSG_NOTICE([Using optimized DP implementation from Jeremy Buhler])
           AC_DEFINE(p7_IMPL_BUHLER)
           ;;

sse)       AC_MSG_NOTICE([Using SSE2 DP implementation for Intel (contrib by Apple)])
           AC_DEFINE(p7_IMPL_SSE)
           ;;

*)         AC_MSG_NOTICE([Ignoring unknown argument to --enable-imp: $impl_choice])
           ;;
esac

CFLAGS="$CFLAGS $impl_cflags"

# For Altivec: test that the compiler is really going to support this.
# (M4 code from Erik Lindahl, Stanford University; 
#  modified for IBM PowerPC xlc tests by Yinhe Cheng, IBM).
#
if test "$impl_choice" = altivec; then
	if test "$ac_cv_build" = powerpc64-unknown-linux-gnu; then
	   cat >> confdefs.h << \_ACEOF
#include "altivec.h"
_ACEOF
        fi
        AC_MSG_CHECKING([whether the compiler supports altivec extensions])
#
# IBM PPC xlc altivec compile test: (Yinhe Cheng, IBM)
#
	if test "$CC" = "xlc"; then  
		AC_TRY_COMPILE([],[
#ifndef __IBM_VACPP_altivec_h
choke_me
#endif
vector signed char X;
	      ],[
	      AC_MSG_RESULT([yes])
	      AC_MSG_NOTICE([Altivec optimization successfully enabled.])
	      ],[
	      AC_MSG_RESULT([no])
	      AC_MSG_ERROR([altivec not supported on this platform]
[]
[(xlc/altivec): Either you aren't on a platform that supports the Altivec]
[instruction set, or your xlc compiler doesn't support Altivec. Rerun configure]
[without --enable-altivec, or update your compiler.])]
)
#
# else, Apple OS/X altivec compile test: (Erik Lindahl, Stanford)
#
	else
	      AC_TRY_COMPILE([],[
#ifndef __VEC__
choke_me
#endif
vector signed char X;
              ],[
              AC_MSG_RESULT([yes])
              AC_MSG_NOTICE([Altivec optimization successfully enabled.])
              ],[
              AC_MSG_RESULT([no])
              AC_MSG_ERROR([altivec not supported on this platform]
[]
[Either you aren't on a PowerPC platform that supports Altivec instructions, or your]
[compiler doesn't support Altivec. Rerun configure without the --enable-altivec]
[flag, or update to a recent compiler like gcc-3.1 (included in the latest Apple]
[developer tools).])]
        )                       
        fi
fi


# --enable-mpi       Enable MPI parallelization
# 
# Sets MPILIBS output variable
# Sets CC to appropriate compiler wrapper for MPI (mpicc, usually)
# Defines HAVE_MPI preprocessor variable.
#
AC_ARG_ENABLE(mpi, 
[AS_HELP_STRING([--enable-mpi], [enable MPI parallelization])],
[case $enable_mpi in
   yes) ACX_MPI([
           AC_MSG_NOTICE([enabled optional MPI parallelization])
	   AC_DEFINE(HAVE_MPI)
	   CC=$MPICC
          ], [
           AC_MSG_ERROR([MPI library not found])
          ])
	 ;;
   no)   AC_MSG_NOTICE([MPI parallelization disabled])
         ;;
   *)    echo "Ignoring unknown argument to --enable-mpi: $enable_mpi"
	 ;;
esac])


# --enable-threads      Enable POSIX multithreading
#
# Uses ACX_PTHREAD macro from the GNU macro archive.
# Back to my code to finish enabling pthreads...
# Define these preprocessor variables:
#     HMMER_THREADS
#     HAVE_PTHREAD_SETCONCURRENCY
#     HAVE_PTHREAD_ATTR_SETSCOPE
#
AC_ARG_ENABLE(threads, 
[  --enable-threads        enable POSIX multithreading support],
[case $enable_threads in 
   yes) AC_MSG_NOTICE([enabled POSIX multithreading support])
        AC_CHECK_FUNCS(pthread_setconcurrency)
        AC_CHECK_FUNCS(pthread_attr_setscope)
	ACX_PTHREAD(AC_DEFINE(HMMER_THREADS))
	;;
   no)  AC_MSG_NOTICE([POSIX multithreading support disabled])
        ;;
   *)   echo "Ignoring unknown argument to --disable-threads: $enable_threads"
	;;
esac])

# --with-gsl        - enable hooks into the GSL (GNU Scientific Library)
#
# The (action-if-given) is blank: using --with-gsl sets $with_gsl to 'yes',
# and --without-gsl sets $with_gsl to 'no'. If neither is given,
# $with_gsl is set to 'check', and we'll try to use GSL anyway if we
# can find it.
# 
# We interpret the $with_gsl setting later, in the "checks for libraries" section.
AC_ARG_WITH([gsl],
            [AS_HELP_STRING([--with-gsl],
                           [use the GSL, GNU Scientific Library (default is no)])],
            [if test "x$withval" != xno ; then 
               AC_MSG_NOTICE([Enabling hooks into the GNU Scientific Library (GSL)])
             fi
            ],
	    [with_gsl=no])


# Checks for programs
#

# Set output variable "RANLIB"
AC_PROG_RANLIB

# Set output variable "AR"
AC_PATH_PROG([AR], [ar], [:], [$PATH:/usr/ccs/bin:/usr/xpg4/bin])

# Set the output variable EXEC_DEPENDENCY according
# to whether we're using GNU make or a SYSV make.
CHECK_GNU_MAKE



# Checks for libraries
#
LIBGSL=
AS_IF([test "x$with_gsl" != xno],
      [AC_CHECK_LIB([gsl], [gsl_expm1], 
           [AC_SUBST([LIBGSL], ["-lgsl -lgslcblas"])
            AC_DEFINE([HAVE_LIBGSL], [1], [Define if you have libgsl])
           ],
           [if test "x$with_gsl" != xcheck; then
             AC_MSG_FAILURE(
               [--with-gsl was given, but GSL library was not found])
            fi
           ],
           [-lgslcblas]
        )])

# Checks for headers
#
AC_CHECK_HEADERS([sys/types.h])
AC_CHECK_HEADERS([unistd.h])
AC_CHECK_HEADERS([stdint.h])
AC_CHECK_HEADERS([inttypes.h])

# Checks for types
#
AC_C_BIGENDIAN
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_OFF_T

# Checks for functions, defining HAVE_FOO when foo is found
#
AC_CHECK_FUNCS(mkstemp)
AC_CHECK_FUNCS(popen)
AC_CHECK_FUNCS(strcasecmp)
AC_CHECK_FUNCS(times)
AC_FUNC_FSEEKO
AC_CHECK_FUNCS(ntohs, , AC_CHECK_LIB(socket, ntohs))
AC_CHECK_FUNCS(ntohl, , AC_CHECK_LIB(socket, ntohl))
AC_CHECK_FUNCS(htons, , AC_CHECK_LIB(socket, htons))
AC_CHECK_FUNCS(htonl, , AC_CHECK_LIB(socket, htonl))

#
# 11. System services
#
AC_SYS_LARGEFILE

# Done.
# Config subdirs and files (except main Makefile, which we defer)
#

AC_CONFIG_HEADERS(src/p7_config.h)
AC_CONFIG_FILES(src/Makefile)
AC_CONFIG_FILES(h3/Makefile)
AC_CONFIG_HEADERS(h3/p7_config.h)
AC_CONFIG_FILES(testsuite/Makefile)

AC_CONFIG_HEADERS([easel/esl_config.h])
AC_CONFIG_FILES([easel/miniapps/Makefile])
AC_CONFIG_FILES([easel/testsuite/Makefile])
AC_CONFIG_FILES([easel/Makefile])


################################################################
# 13. AC_OUTPUT
################################################################
# remember that sneaky if statement that wraps the source
# configuration? well, it ends here.
else
  AC_MSG_NOTICE([This is a precompiled binary distribution - using abbreviated config])
fi
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

echo "

HMMER configuration:

     DP implementation:    ${impl_choice}
"




