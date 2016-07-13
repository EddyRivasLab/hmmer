# Process this file with autoconf to produce the HMMER configure script.
#
# This is synchronized with Easel's configure script. HMMER's
# configure script configures Easel in addition to itself. This saves
# having to do a separate ./configure in Easel, but also means that
# everything that appears in Easel's configure script must also appear
# here.
#
# Some reminders to save re-reading autoconf manual for the n'th time:
#    output variables:
#      -  are defined here as normal shell variables, e.g. FOO="my string"
#      -  are made into output variables by calling AC_SUBST(FOO)
#      -  any occurrence of @FOO@ in an output file is then substituted
#         This only happens in files we assign w/ AC_CONFIG_FILES;
#         such as our Makefile.in's, for example.
#
#   C preprocessor symbols:
#      -  are defined here by calling AC_DEFINE(FOO) or AC_DEFINE(FOO, [42])
#      -  then #undef FOO in easel.h.in becomes #define FOO, #define FOO 42
#         This only happens in header files that we assign
#         w/ AC_CONFIG_HEADERS, such as p7_config.h.in
#
# Contents:   
#    1. autoconf requirements
#    2. AC_INIT
#    3. Info on the package
#    4. Process ./configure commandline options
#    5. Checks for programs, including ${CC}, ${CFLAGS}
#    6. Checks for libraries, such as GSL
#    7. Checks for headers, including SSE headers
#    8. Checks for types.
#    9. Checks for structures.
#   10. Checks for compiler characteristics
#   11. Checks for functions.
#   12. Checks for system services
#   13. AC_CONFIG_FILES
#   14. AC_OUTPUT
#   15. Copyright and license information
#
# Order obeys autoconf manual, "standard configure.ac layout".



################################################################
# 1. autoconf requirements
################################################################
# autoconf 2.61 has a bug in AC_FUNC_FSEEKO; make sure we don't use it.
AC_PREREQ(2.63)


################################################################
# 2. AC_INIT
################################################################

AC_INIT(HMMER, 4.0dev, sean@eddylab.org, hmmer)
AC_MSG_NOTICE([Configuring HMMER for your system.])

# remember if the user is overriding CFLAGS
if test x"$CFLAGS" = x; then 
  sre_cflags_env_set=no
else
  sre_cflags_env_set=yes
fi

################################################################
# 3. Info on the package
################################################################
#
# The four AC_INIT args above set these output variables and preprocessor symbols:
#     PACKAGE_NAME      <package>     e.g. "HMMER"
#     PACKAGE_VERSION   <version>     e.g. "4.0dev"
#     PACKAGE_BUGREPORT <bug-report>  e.g. "sean@eddylab.org"
#     PACKAGE_TARNAME   <tarname>     e.g. "hmmer"
# From them, it derives one more:
#     PACKAGE_STRING    <package> <version>, e.g. "HMMER 4.0dev"
# We also define additional variables:
#     HMMER_DATE        release date: e.g. "January 2007"
#     HMMER_COPYRIGHT   one-line copyright string
#     HMMER_LICENSE     one-line license string
#     HMMER_LICENSETAG  which license to bundle from Licenses/ subdirectory.
#     HMMER_VERSION     copy of version code, e.g. "3.0"
#     HMMER_URL         URL home for HMMER:   e.g. "hmmer.org"
# And we have to define the relevant package variables for Easel as well.
#
# We avoid using AC_INIT's PACKAGE_ variables anywhere, because we want to be able
# to use HMMER as a library inside other packages, with no name clashes.
################################################################

HMMER_DATE="Apr 2015"
HMMER_COPYRIGHT="Copyright (C) 2015 Howard Hughes Medical Institute."
HMMER_LICENSE="Freely distributed under the open source BSD license."
HMMER_VERSION=$PACKAGE_VERSION
HMMER_URL="hmmer.org"
HMMER_ESLDIR="lib/easel"
HMMER_SADIR="lib/libdivsufsort"

EASEL_DATE="Apr 2015"
EASEL_COPYRIGHT="Copyright (C) 2015 Howard Hughes Medical Institute."
EASEL_LICENSE="Freely distributed under the open source BSD license."
EASEL_VERSION="h$PACKAGE_VERSION"
EASEL_URL="http://bioeasel.org/"

# Output variables  (AC_OUTPUT replaces @var@ in input files, such as Makefiles)
AC_SUBST(HMMER_DATE)
AC_SUBST(HMMER_COPYRIGHT)
AC_SUBST(HMMER_LICENSE)
AC_SUBST(HMMER_VERSION)
AC_SUBST(HMMER_URL)
AC_SUBST(HMMER_ESLDIR)
AC_SUBST(HMMER_SADIR)

AC_SUBST(EASEL_DATE)
AC_SUBST(EASEL_COPYRIGHT)
AC_SUBST(EASEL_LICENSE)
AC_SUBST(EASEL_VERSION)
AC_SUBST(EASEL_URL)

# Preprocessor symbols (replace #undefs in p7_config.h)
AC_DEFINE_UNQUOTED(HMMER_DATE,      "$HMMER_DATE")
AC_DEFINE_UNQUOTED(HMMER_COPYRIGHT, "$HMMER_COPYRIGHT")
AC_DEFINE_UNQUOTED(HMMER_LICENSE,   "$HMMER_LICENSE")
AC_DEFINE_UNQUOTED(HMMER_VERSION,   "$HMMER_VERSION")
AC_DEFINE_UNQUOTED(HMMER_URL,       "$HMMER_URL")

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
#if ((! test -d binaries) && test -d src); then 
AC_MSG_NOTICE([Full HMMER source distribution - using full configuration])



################################################################
# 4. Process the ./configure command line
################################################################

# --enable-debugging      - set basic debugging (level 1)
# --enable-debugging=x    - set debugging level to <x> (1-3)
#
# At all levels, including 0, replaces CFLAGS w/ "-g -Wall" (so it assumes gcc).
# Sets p7_DEBUGGING preprocessor symbol, which compiles in debugging support.
# Also sets eslDEBUGLEVEL level, 0-3.
#
AC_ARG_ENABLE(debugging,
[AS_HELP_STRING([--enable-debugging],[include debugging code])
AS_HELP_STRING([--enable-debugging=x],[also set diagnostics verbosity level to <x> (1-3)])],
enable_debugging=$enableval, enable_debugging=no)

case $enable_debugging in
   yes)  AC_DEFINE(eslDEBUGLEVEL, 1);;
     1)  AC_DEFINE(eslDEBUGLEVEL, 1);;
     2)  AC_DEFINE(eslDEBUGLEVEL, 2);;
     3)  AC_DEFINE(eslDEBUGLEVEL, 3);;
    no)  AC_DEFINE(eslDEBUGLEVEL, 0);;
     *)  AC_MSG_ERROR([Unknown argument to --enable-debugging: $enable_debugging]);;
esac

# --enable-exact-logsum  - turn on p7_LOGSUM_SLOWEXACT flag (should only be used in debugging)
# 
AC_ARG_ENABLE(exact-logsum,
  [AS_HELP_STRING( [--enable-exact-logsum],[use slow/exact lse2 log-sum-exp calculations])],
  ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  AC_MSG_NOTICE([Compiling with slow/exact lse2 log-sum-exp calculations.])
  AC_DEFINE(p7_LOGSUM_SLOWEXACT)
fi


# --enable-gcov      - compile for code coverage testing
#
# Replaces CC with "gcc" and CFLAGS with "-g -Wall -fprofile-arcs -ftest-coverage".
# Running programs (notably the test suite) will then generate .gcda files,
# which the gcov tool can read.
#
AC_ARG_ENABLE(gcov,[AS_HELP_STRING([--enable-gcov],[compile for code coverage testing])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  AC_MSG_NOTICE([Compiling with gcov instrumentation.])
  CC="gcc"
  CFLAGS="-g -Wall -fprofile-arcs -ftest-coverage"
  sre_cflags_env_set=yes
fi


# --enable-gprof     - compile for gcc code profiling
#
# Replaces CC with "gcc" and CFLAGS with "-O -g -pg"
# Running programs generates a profile data file "gmon.out".
#
AC_ARG_ENABLE(gprof,[AS_HELP_STRING([--enable-gprof],[compile for gcc code profiing])], ok=$enableval, ok=no)
if test "$ok" = "yes"; then
  AC_MSG_NOTICE([Compiling with gprof instrumentation.])
  CC="gcc"
  CFLAGS="-O -g -pg"
  sre_cflags_env_set=yes
fi


# --enable-mpi       Enable MPI parallelization
AC_ARG_ENABLE(mpi,[AS_HELP_STRING([--enable-mpi],[enable MPI parallelization])], enable_mpi=$enableval, enable_mpi=no)

# --enable-xlc-arch  Set the code generation for a particular architecture
#
# The xlc compiler on AIX seems to need the -qarch flag to be set for a processor
# that is VMX enabled.  Otherwise a warning is issued that states the -qenablevmx
# and -qaltivec flags are ignored.
AX_COMPILER_VENDOR
if test "${ax_cv_c_compiler_vendor}" = "ibm"; then
  AC_ARG_WITH(xlc-arch, 
              [AC_HELP_STRING([--with-xlc-arch=<arch>], [specify architecture <arch> for xlc -qarch])],
              XLC_ARCH=$withval, 
              XLC_ARCH=no)
fi

# --enable-threads      Enable POSIX multithreading
#
# Uses ACX_PTHREAD macro from the GNU macro archive.
# Back to my code to finish enabling pthreads...
# Define these preprocessor variables:
#     HMMER_THREADS
#     HAVE_PTHREAD_SETCONCURRENCY
#     HAVE_PTHREAD_ATTR_SETSCOPE
#
AC_ARG_ENABLE([threads],
              [AS_HELP_STRING([--enable-threads],
                 [enable POSIX multithreading support (default is check)])],
	      [],
              [enable_threads=check])

case $enable_threads in 
   yes)   check_threads=yes
          ;;
   no)    check_threads=no
          ;;
   check) check_threads=yes
          ;;
   *)     check_threads=no
          AC_MSG_NOTICE([Ignoring unknown argument to --enable-threads: $enable_threads])
          ;;
esac
 
AS_IF([test "x$check_threads" != xno],
      [ACX_PTHREAD([AC_DEFINE(HMMER_THREADS)
                    AC_DEFINE(HAVE_PTHREAD)
                    AC_SUBST(PTHREAD_LIBS)
                    AC_SUBST(PTHREAD_CFLAGS)],
                   [if test "x$enable_threads" != xcheck; then
                      AC_MSG_FAILURE([POSIX multithreading failed])
                    fi])])

#AC_ARG_ENABLE(threads, 
#[  --enable-threads        enable POSIX multithreading support],
#[case $enable_threads in 
#   yes) AC_MSG_NOTICE([enabled POSIX multithreading support])
#        AC_CHECK_FUNCS(pthread_setconcurrency)
#        AC_CHECK_FUNCS(pthread_attr_setscope)
#	ACX_PTHREAD(AC_DEFINE(HMMER_THREADS))
#	;;
#   no)  AC_MSG_NOTICE([POSIX multithreading support disabled])
#        ;;
#   *)   echo "Ignoring unknown argument to --disable-threads: $enable_threads"
#	;;
#esac])

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

# Enable NEON extension on ARM processors
AC_ARG_ENABLE(neon,[AS_HELP_STRING([--enable-neon],[enable NEON extension])],          enable_neon=$enableval,   enable_neon=check)
################################################################
# 5. Checks for programs, including ${CC} and its ${CFLAGS}.
#################################################################
AC_PROG_CC
AC_PROG_INSTALL
# MPI :  set @CC@ to mpicc;
#        set @MPILIBS@ if needed (usually not; mpicc deals w/ it);
#        defines HAVE_MPI. 
if test "$enable_mpi" = "yes"; then
   ACX_MPI([], AC_MSG_ERROR([MPI library not found for --enable-mpi]))
   CC=$MPICC
fi
AC_PROG_CPP

# AC_PROG_CC_STDC wants to append -c99 to icc; icc doesn't like this.
# To autodetect icc, we use the AX_COMPILER_VENDOR macro.
AX_COMPILER_VENDOR
if test "$ax_cv_c_compiler_vendor" != "intel"; then
  AC_PROG_CC_STDC
fi
AC_PROG_RANLIB
AC_PATH_PROG([AR], [ar], [:], [$PATH:/usr/ccs/bin:/usr/xpg4/bin])


# CHECK_GNU_MAKE sets EXEC_DEPENDENCY for SYSV vs. GNU make,
# to support Makefile lines of ${PROGS}: ${EXEC_DEPENDENCY}
CHECK_GNU_MAKE


# If the base $CFLAGS isn't already set, select what it should be.
# We've used AX_CC_MAXOPT in the past (from FFTW), but it's
# superaggressive and problematic in many respects. Finally gave up on
# it. -O3 is enough.
if test "$enable_debugging" != "no"; then
   AC_DEFINE(p7_DEBUGGING,  1)
   if test "$GCC" = "yes" && test "$sre_cflags_env_set" = "no"; then
      CFLAGS="-g -Wall"
   fi
else
   if test "$sre_cflags_env_set" = "no"; then
      CFLAGS="-O3"	
   fi
fi

# Figure out how to generate PIC code for shared libraries
# This sets @PIC_FLAGS@ output variable for Easel.
ESL_PIC_FLAGS

# ARM and NEON SIMD flags
if test "$enable_neon" = "yes" || test "$enable_neon" = "check"; then
   if test "x$SIMD_CFLAGS" = x; then
     case $ax_cv_c_compiler_vendor in
       gnu)   AX_CHECK_COMPILER_FLAGS(-mfpu=neon,[SIMD_CFLAGS="-mfpu=neon"],[]);;
       *)     ;;
     esac
   fi
fi
if test "$enable_neon" = "check"; then
   sre_save_CFLAGS="$CFLAGS"
   sre_save_CPPFLAGS="$CPPFLAGS"
   CFLAGS="$CFLAGS $SIMD_CFLAGS"
   CPPFLAGS="$CPPFLAGS $SIMD_CFLAGS"
   AC_CHECK_HEADER([arm_neon.h],[
     AC_MSG_CHECKING([if platform can run NEON instructions])
     AC_RUN_IFELSE([AC_LANG_PROGRAM( [[#include <arm_neon.h>]],
                              [[ int32x4_t x = vdupq_n_s32(0);]])],
                   [AC_MSG_RESULT([yes])
                     enable_neon=yes],
                   [AC_MSG_RESULT([no])])])
   CFLAGS="$sre_save_CFLAGS"
   CPPFLAGS="$sre_save_CPPFLAGS"
fi

# Determine target if targeting is specified.
# Do this by attempting to execute code with different 
# cominations of NEON/VFP compiler options. This allows 
# us to determine the best compile options for each 
# ARM architecture. 
if test "$enable_neon" = "yes" || test "$enable_neon" = "check"; then
  cortex_a5="-mcpu=cortex-a5 -mfpu=neon-fp16"
  cortex_a7="-mcpu=cortex-a7 -mfpu=neon-vfpv4"
  cortex_a8="-mcpu=cortex-a8 -mfpu=neon"
  cortex_a9="-mcpu=cortex-a9 -mfpu=neon-fp16"
  cortex_a15="-mcpu=cortex-a15 -mfpu=neon-vfpv4"
  sre_save_CFLAGS="$CFLAGS"
  sre_save_CPPFLAGS="$CPPFLAGS"
  CFLAGS="$CFLAGS $cortex_a5"
  CPPFLAGS="$CPPFLAGS $cortex_a5"     
  AC_RUN_IFELSE([AC_LANG_PROGRAM( [[]],[[int a = 32;]])],
                    [vector_flags="$cortex_a5"], [])
  CFLAGS="$sre_save_CFLAGS $cortex_a7"  
  CPPFLAGS="$sre_save_CPPFLAGS $cortex_a7"
  AC_RUN_IFELSE([AC_LANG_PROGRAM( [[]],[[int a = 32;]])],
                    [vector_flags="$cortex_a7"], [])
  CFLAGS="$sre_save_CFLAGS $cortex_a8"  
  CPPFLAGS="$sre_save_CPPFLAGS $cortex_a8"
  AC_RUN_IFELSE([AC_LANG_PROGRAM( [[]],[[int a = 32;]])],
                    [vector_flags="$cortex_a8"], [])
  CFLAGS="$sre_save_CFLAGS $cortex_a9"  
  CPPFLAGS="$sre_save_CPPFLAGS $cortex_a9"
  AC_RUN_IFELSE([AC_LANG_PROGRAM( [[]],[[int a = 32;]])],
                    [vector_flags="$cortex_a9"], [])
  CFLAGS="$sre_save_CFLAGS $cortex_a15"  
  CPPFLAGS="$sre_save_CPPFLAGS $cortex_a15"
  AC_RUN_IFELSE([AC_LANG_PROGRAM( [[]],[[int a = 32;]])],
                    [vector_flags="$cortex_a15"], [])
  
  CFLAGS="$sre_save_CFLAGS"
  CPPFLAGS="$sre_save_CPPFLAGS"

fi
if test "x$vector_flags" = x; then
  SIMD_CFLAGS="-mfpu=neon" 
fi


if test "y$vector_flags" != y; then
  SIMD_CFLAGS="$vector_flags"
  enable_neon=yes
fi

if test "$enable_neon" = yes; then
  AC_DEFINE(HAVE_NEON,1,[Support NEON instructions])
fi

AC_SUBST(SIMD_CFLAGS)


# Define HAVE_GZIP if gzip is in $PATH (or if HAVE_GZIP is already set)
AC_PATH_PROG(HAVE_GZIP, "gzip", "no")
if test "${HAVE_GZIP}" = "no"; then
  AC_MSG_WARN([gzip not found])
else
  AC_DEFINE(HAVE_GZIP)
fi



################################################################
# 6. Checks for libraries
#################################################################
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

# Easel stopwatch high-res timer may try to use clock_gettime,
# which may be in librt
AC_SEARCH_LIBS(clock_gettime, [rt posix4])


################################################################
# 7. Checks for headers
################################################################
# Defines HAVE_SYS_TYPES_H, HAVE_STDINT_H, etc.
AC_CHECK_HEADERS([ \
  endian.h\
  inttypes.h\
  stdint.h\
  unistd.h\
  sys/types.h\
  netinet/in.h
])  

# Check for sysctl.h separately.  On OpenBSD, it requires
# <sys/param.h> and autoconf needs special logic to deal w. this as
# follows.
AC_CHECK_HEADERS([sys/param.h])
AC_CHECK_HEADERS([sys/sysctl.h], [], [], 
[[#ifdef HAVE_SYS_PARAM_H
#include <sys/param.h>
#endif
]])


# Vector-specific headers. 
AC_CHECK_HEADER([arm_neon.h])

################################################################
# 8. Checks for types
################################################################
#    - Define WORDS_BIGENDIAN on bigendian platforms.
#    - Make sure we have C99 exact-size integer types; 
#      ssi uses 16, 32, and 64-bit ints, and we
#      use 8-bit unsigned chars for digitized sequence.
#    - Make sure we have off_t.
AC_C_BIGENDIAN
AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_OFF_T

################################################################
# 9. Checks for structures - currently none
#################################################################



################################################################
# 10. Checks for compiler characteristics
#################################################################

# __attribute__() tags on function declarations
# HAVE_FUNC_ATTRIBUTE_NORETURN 
#
#   The clang static analyzer can't figure out that some of our
#   varargs-dependent fatal error handlers (esl_fatal(), for example)
#   cannot return. To tell it so, we take advantage of __attribute__
#   tags on function declarations, a non-ISO gcc extension, when
#   available. gcc, clang, and other gcc-like compilers support this.
# 
AX_GCC_FUNC_ATTRIBUTE(noreturn)


################################################################
# 11. Checks for functions, defining HAVE_FOO when foo is found
#################################################################
AC_CHECK_FUNCS(mkstemp)
AC_CHECK_FUNCS(popen)
AC_CHECK_FUNCS(strcasecmp)
AC_CHECK_FUNCS(times)
AC_CHECK_FUNCS(getpid)
AC_CHECK_FUNCS(sysctl)
AC_CHECK_FUNCS(sysconf)
AC_CHECK_FUNCS(getcwd)
AC_CHECK_FUNCS(stat)
AC_CHECK_FUNCS(fstat)
AC_CHECK_FUNCS(erfc)

AC_CHECK_FUNCS(ntohs, , AC_CHECK_LIB(socket, ntohs))
AC_CHECK_FUNCS(ntohl, , AC_CHECK_LIB(socket, ntohl))
AC_CHECK_FUNCS(htons, , AC_CHECK_LIB(socket, htons))
AC_CHECK_FUNCS(htonl, , AC_CHECK_LIB(socket, htonl))
AC_FUNC_FSEEKO

#################################################################
# 12. System services
#################################################################
AC_SYS_LARGEFILE


################################################################
# 13. Config subdirs and files (except main Makefile, which we defer)
#################################################################
if test -d documentation/release-notes; then
  AC_CONFIG_FILES(documentation/release-notes/LICENSE.sh)
fi

AC_CONFIG_FILES([
  documentation/Makefile
  documentation/man/Makefile
  documentation/userguide/Makefile
  src/Makefile
  src/Makefile-subdirs.mk
  src/base/Makefile
  src/build/Makefile
  src/daemon/Makefile
  src/dp_reference/Makefile
  src/dp_sparse/Makefile
  src/dp_vector/Makefile
  src/experiments/Makefile
  src/misc/Makefile
  src/programs/Makefile
  src/sandbox/Makefile
  src/search/Makefile
  src/utilities/Makefile
  testsuite/Makefile
  benchmarks/Makefile
  benchmarks/profmark/Makefile
  benchmarks/speed/Makefile
  lib/easel/LICENSE.sh
  lib/easel/miniapps/Makefile
  lib/easel/testsuite/Makefile
  lib/easel/documentation/Makefile
  lib/easel/Makefile
  lib/libdivsufsort/Makefile
])
  
AC_CONFIG_HEADERS([
  src/p7_config.h
  lib/easel/esl_config.h
  lib/libdivsufsort/divsufsort.h
])

################################################################
# 14. et voila!
################################################################
# remember that sneaky if statement that wraps the source
# configuration? well, it ends here.
#else
#  AC_MSG_NOTICE([This is a precompiled binary distribution - using abbreviated config])
#fi
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

echo "

HMMER configuration:
     compiler:             ${CC} ${CFLAGS} ${SIMD_CFLAGS} ${PTHREAD_CFLAGS} ${PIC_FLAGS}
     host:                 $host
     linker options:       ${LDFLAGS}
     libraries:            ${LIBGSL} ${MPILIBS} ${PTHREAD_LIBS} ${LIBS}
"


#################################################################
# @LICENSE@
# 
# SVN $Id$
# SVN $URL$ 
#################################################################

