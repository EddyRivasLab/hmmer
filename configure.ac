# Process this file with autoconf to produce the HMMER configure script.
#
# Note that this is synchronized with Easel's configure script,
# so that HMMER configures Easel in addition to itself; this saves
# having to do a separate ./configure in Easel. That is,
# everything that appears in Easel's configure script must also
# appear here. 
#
# reminders to save re-reading autoconf manual for the n'th time:
#    output variables:
#      -  are defined here as normal shell variables, e.g. FOO="my string"
#      -  are made into output variables by calling AC_SUBST(FOO)
#      -  any occurrence of @FOO@ in an output file is then substituted
#         I think this only happens in files we assign w/ AC_CONFIG_FILES;
#         that's the Makefile.in's.
#
#   C preprocessor symbols:
#      -  are defined here by calling AC_DEFINE(FOO) or AC_DEFINE(FOO, [42])
#      -  then #undef FOO in easel.h.in becomes #define FOO, #define FOO 42
#         I think this only happens in header files that we assign
#         w/ AC_CONFIG_HEADERS -- which means, easel.h.in
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


# Autoconf 2.61 has a bug in AC_FUNC_FSEEKO; don't use it.
AC_PREREQ(2.63)
AC_INIT(HMMER, 3.0b2, eddys@janelia.hhmi.org, hmmer)
AC_MSG_NOTICE([Configuring HMMER for your system.])

# remember if the user is overriding CFLAGS
if test x"$CFLAGS" = x; then 
  sre_cflags_env_set=no
fi

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

HMMER_DATE="June 2009"
HMMER_COPYRIGHT="Copyright (C) 2009 Howard Hughes Medical Institute."
HMMER_LICENSE="Freely distributed under the GNU General Public License (GPLv3)."
HMMER_VERSION=$PACKAGE_VERSION
HMMER_URL="http://hmmer.org/"

EASEL_DATE="June 2009"
EASEL_COPYRIGHT="Copyright (C) 2009 Howard Hughes Medical Institute."
EASEL_LICENSE="Freely distributed under the Janelia Farm Software License."
EASEL_VERSION="h$PACKAGE_VERSION"
EASEL_URL="http://bioeasel.org/"

# Output variables  (AC_OUTPUT replaces @var@ in input files, such as Makefiles)
AC_SUBST(HMMER_DATE)
AC_SUBST(HMMER_COPYRIGHT)
AC_SUBST(HMMER_LICENSE)
AC_SUBST(HMMER_VERSION)
AC_SUBST(HMMER_URL)

AC_SUBST(EASEL_DATE)
AC_SUBST(EASEL_COPYRIGHT)
AC_SUBST(EASEL_LICENSE)
AC_SUBST(EASEL_VERSION)
AC_SUBST(EASEL_URL)

# Preprocessor symbols (replace #undefs in p7config.h)
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
# Process the ./configure command line
################################################################

# --enable-debugging      - set basic debugging (level 0)
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
   yes)  AC_DEFINE(eslDEBUGLEVEL, 0);;
     1)  AC_DEFINE(eslDEBUGLEVEL, 1);;
     2)  AC_DEFINE(eslDEBUGLEVEL, 2);;
     3)  AC_DEFINE(eslDEBUGLEVEL, 3);;
    no)  AC_DEFINE(eslDEBUGLEVEL, 0);;
     *)  AC_MSG_ERROR([Unknown argument to --enable-debugging: $enable_debugging]);;
esac


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


# --enable-sse    Enable the SSE implementation.
# --enable-vmx    Enable the Altivec/VMX implementation.
# --enable-dummy  Enable the portable implementation.
#
# These set $impl_choice to one of sse, vmx, dummy, or none.
#
AC_ARG_ENABLE(sse,[AS_HELP_STRING([--enable-sse],[enable SSE implementation])],          enable_sse=$enableval,   enable_sse=no)
AC_ARG_ENABLE(vmx,[AS_HELP_STRING([--enable-vmx],[enable Altivec/VMX implementation])],  enable_vmx=$enableval,   enable_vmx=no)
AC_ARG_ENABLE(dummy,[AS_HELP_STRING([--enable-dummy],[enable portable implementation])], enable_dummy=$enableval, enable_dummy=no)
if test "$enable_sse" = "yes"; then
   if test "$enable_vmx" = "yes" || test "$enable_dummy" = "yes"; then
      AC_MSG_ERROR([Select only one implementation of sse, vmx, or dummy])
   fi
   impl_choice=sse
elif test "$enable_vmx" = "yes"; then
   if test "$enable_sse" = "yes" || test "$enable_dummy" = "yes"; then
      AC_MSG_ERROR([Select only one implementation of sse, vmx, or dummy])
   fi
   impl_choice=vmx
elif test "$enable_dummy" = "yes"; then
   if test "$enable_sse" = "yes" || test "$enable_vmx" = "yes"; then
      AC_MSG_ERROR([Select only one implementation of sse, vmx, or dummy])
   fi
   impl_choice=dummy
else
   impl_choice=none   	
fi


# --enable-mpi       Enable MPI parallelization
AC_ARG_ENABLE(mpi,[AS_HELP_STRING([--enable-mpi],[enable MPI parallelization])], enable_mpi=$enableval, enable_mpi=no)


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
                    AC_DEFINE(HAVE_PTHREADS)
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


# If the user didn't force an implementation on us, select a default
if test "$impl_choice" = "none"; then
  case $host in 
    ia64-*-*)      impl_choice=sse;;
    i?86-*-*)      impl_choice=sse;;
    x86*-*-*)      impl_choice=sse;;
    powerpc*-*-*)  impl_choice=vmx;;
    *)             impl_choice=dummy;;
  esac
fi

################################################################
# Checks for programs, including choice of compiler ${CC} and ${CFLAGS}.
#
AC_PROG_CC
AC_PROG_CPP
AX_COMPILER_VENDOR
# AC_PROG_CC_STDC wants to append -c99 to icc; icc doesn't like this.
if test "$ax_cv_c_compiler_vendor" != "intel"; then
  AC_PROG_CC_STDC
fi
AC_PROG_RANLIB
AC_PATH_PROG([AR], [ar], [:], [$PATH:/usr/ccs/bin:/usr/xpg4/bin])
CHECK_GNU_MAKE

# Select our optimization level in CFLAGS.
if test "$enable_debugging" != "no"; then
   AC_DEFINE(p7_DEBUGGING,  1)
   if test "$GCC" = "yes" && test "$sre_cflags_env_set" = "no"; then
      CFLAGS="-g -Wall"
   fi
else
   sre_save_ac_test_cflags="$ac_test_CFLAGS"
   if test "$sre_cflags_env_set" = "no"; then
     ac_test_CFLAGS=""
   fi
   AX_CC_MAXOPT
   ac_test_CFLAGS="$sre_save_ac_test_cflags"
   # That was a workaround. AX_CC_MAXOPT apparently thinks that
   # $ac_test_CFLAGS only gets set if ${CFLAGS} was set by the
   # user (in the environment), but that's not so; it gets set
   # by AC_PROG_CC. Instead of rewriting AX_CC_MAXOPT, we fool
   # it.
fi

# MPI :  set @CC@ to mpicc, sets @MPILIBS@, defines HAVE_MPI. 
if test "$enable_mpi" = "yes"; then
   ACX_MPI(,AC_MSG_ERROR([MPI library not found for --enable-mpi]))
   CC=$MPICC
   if test "$ax_cv_c_compiler_vendor" != "intel"; then
     AC_PROG_CC_STDC
   fi
fi

# If we're using SSE, figure out our @SIMD_CFLAGS@
if test "$impl_choice" = "sse" && test "x$SIMD_CFLAGS" = x; then
   case $ax_cv_c_compiler_vendor in
   gnu)   AX_CHECK_COMPILER_FLAGS(-msse2,[SIMD_CFLAGS="-msse2"],[impl_choice=dummy]);;
   intel) ;;
   *)     ;;
   esac
fi
if test "$enable_sse" = "yes" && test "$impl_choice" = dummy; then
   AC_MSG_ERROR([SSE instructions are apparently not available; --enable-sse failed])
fi

# If we're using VMX/Altivec, figure out our @SIMD_CFLAGS@
if test "$impl_choice" = "vmx"; then
   if test "x$SIMD_CFLAGS" = x; then
      if test "$GCC" = "yes"; then
       	# -DFAKE__VEC__ is a workaround because gcc-3.3 does not
	# #define __VEC__ with -maltivec.
        AX_CHECK_COMPILER_FLAGS(-faltivec, [SIMD_CFLAGS="-faltivec"],
	    [AX_CHECK_COMPILER_FLAGS(-maltivec -mabi=altivec,
		[SIMD_CFLAGS="-maltivec -mabi=altivec"],
#		[SIMD_CFLAGS="-maltivec -mabi=altivec -DFAKE__VEC__"],
	        [AX_CHECK_COMPILER_FLAGS(-fvec, [SIMD_CFLAGS="-fvec"],
		   [AC_MSG_ERROR([Need a version of gcc with -maltivec])])])])
      elif test "${ax_cv_c_compiler_vendor}" = "ibm"; then
        AX_CHECK_COMPILER_FLAGS(-qaltivec, [SIMD_CFLAGS="-qaltivec"],[])
      else
        AX_CHECK_COMPILER_FLAGS(-faltivec, [SIMD_CFLAGS="-faltivec"],[])
      fi
   fi
fi

# check if the SSE2 implementation support cast functions
if test "$impl_choice" = "sse"; then
  AC_MSG_CHECKING([checking compiler support for sse2 cast functions])
  ssecast=no
  AC_LANG(C)
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <emmintrin.h>]],
                                     [[__m128 a;
__m128i b;
b = _mm_castps_si128(a);
a = _mm_castsi128_ps(b);]])],
                    [AC_DEFINE([HAVE_SSE2_CAST])
                     ssecast=yes])
  AC_MSG_RESULT([$ssecast])
fi

# Now, we can enable the appropriate optimized implementation.
case "$impl_choice" in 
sse)  AC_MSG_NOTICE([Activating Intel/AMD SSE optimized DP implementation])
      AC_DEFINE([HAVE_SSE2])
      AC_DEFINE([p7_IMPL_SSE])
      AC_SUBST([HMMERIMPLLIB], ["impl_sse/libhmmerimpl.a"])
      IMPL_CHOICE="sse"
      ;;

vmx)  AC_MSG_NOTICE([Activating Altivec/VMX optimized DP implementation])
      AC_DEFINE([HAVE_VMX])
      AC_DEFINE([p7_IMPL_VMX])
      AC_SUBST([HMMERIMPLLIB], ["impl_vmx/libhmmerimpl.a"])
      IMPL_CHOICE="vmx"
      ;;

*)    AC_MSG_NOTICE([No optimized DP implementation available])
      AC_DEFINE([p7_IMPL_DUMMY])
      IMPL_CHOICE="dummy"
      ;;
esac
AC_SUBST(IMPL_CHOICE)
AC_SUBST(SIMD_CFLAGS)


# Check if the linker supports library groups for recursive libraries
AS_IF([test "x$impl_choice" != xno],
      [AC_MSG_CHECKING([compiler support --start-group])
       LDFLAGS_save=$LDFLAGS
       LDFLAGS="-Wl,--start-group -Wl,--end-group $LDFLAGS"
       AC_LINK_IFELSE(AC_LANG_PROGRAM([],[]),
           [AC_MSG_RESULT([yes])
            AC_SUBST([GROUPHMMERLIBS], ["-Wl,--start-group -lhmmer -lhmmerimpl -Wl,--end-group"])],
           [AC_MSG_RESULT([no])
            AC_SUBST([GROUPHMMERLIBS], ["-lhmmer -lhmmerimpl"])])
       LDFLAGS=$LDFLAGS_save],
      [AC_SUBST([GROUPHMMERLIBS], ["-lhmmer"])])

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
AC_CHECK_HEADERS([sys/types.h sys/sysctl.h unistd.h stdint.h inttypes.h])

# altivec.h requires the simd cflags
# For reasons I don't understand, this needs to come after any other CHECK_HEADERS().
if test "$impl_choice" = "vmx"; then
   sre_save_CFLAGS="$CFLAGS"
   sre_save_CPPFLAGS="$CPPFLAGS"
   CFLAGS="$CFLAGS $SIMD_CFLAGS"
   CPPFLAGS="$CPPFLAGS $SIMD_CFLAGS"
   AC_CHECK_HEADERS([altivec.h])
   CFLAGS="$sre_save_CFLAGS"
   CPPFLAGS="$sre_save_CPPFLAGS"
fi

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
AC_CHECK_FUNCS(getpid)
AC_CHECK_FUNCS(sysctl)
AC_CHECK_FUNCS(sysconf)

AC_CHECK_FUNCS(ntohs, , AC_CHECK_LIB(socket, ntohs))
AC_CHECK_FUNCS(ntohl, , AC_CHECK_LIB(socket, ntohl))
AC_CHECK_FUNCS(htons, , AC_CHECK_LIB(socket, htons))
AC_CHECK_FUNCS(htonl, , AC_CHECK_LIB(socket, htonl))
AC_FUNC_FSEEKO

#
# 11. System services
#
AC_SYS_LARGEFILE


# Done.
# Config subdirs and files (except main Makefile, which we defer)
#
if test -d release-notes; then
  AC_CONFIG_FILES(release-notes/LICENSE.csh)
fi
if test -d documentation/alphaguide; then  
  AC_CONFIG_FILES([documentation/alphaguide/Makefile])
fi
AC_CONFIG_FILES([src/Makefile testsuite/Makefile profmark/Makefile])

AC_CONFIG_FILES([src/impl_${impl_choice}/Makefile])
AC_CONFIG_HEADERS(src/p7_config.h)

AC_CONFIG_HEADERS([easel/esl_config.h])
AC_CONFIG_FILES([easel/miniapps/Makefile])
AC_CONFIG_FILES([easel/testsuite/Makefile])
AC_CONFIG_FILES([easel/Makefile])


################################################################
# 13. AC_OUTPUT
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

     compiler:             ${CC} ${CFLAGS} ${SIMD_CFLAGS} ${PTHREAD_CFLAGS}
     libraries:            ${LIBS} ${LIBGSL} ${PTHREAD_LIBS}
     DP implementation:    ${impl_choice}
"




