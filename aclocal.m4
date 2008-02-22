# aclocal.m4 contains custom macros used for creating HMMER's
# configuration script.
#
#   1. ACX_MPI     Configure MPI 
#
#
# SRE, Sun Apr 22 09:26:38 2007 [Janelia]
# SVN $Id$

#################################################################
# Macro: CHECK_GNU_MAKE
# Usage: CHECK_GNU_MAKE
# Author: John Darrington <j.darrington@elvis.murdoch.edu.au> 
# Modified from the original.
# 
# Sets the format of makefile dependency lines for executables.
#
# We need this because GNU make and SYSV make use different systems
# specifying variables for dependencies: $$@ in sysv, %: %.o in GNU.
# Would love to hear a better way of doing this.
# 
# I use two different conventions in my Makefiles. Sometimes 
# executable "foo" has a file "foo.c" - this is the HMMER, Easel, Infernal convention.
# Sometimes executable "foo" has a file "foo_main.c" - this is
# the SQUID convention. The configure script sets the
# EXEC_DEPENDENCY appropriately: here, HMMER style.
#
# Sets an output variable EXEC_DEPENDENCY. 
# This is used in the src/Makefile.in.
#
AC_DEFUN(CHECK_GNU_MAKE,[ 
  AC_MSG_CHECKING(whether your make is GNU make)
  foundGNUmake='nope, assuming sysv make.' ;
  EXEC_DEPENDENCY=[\$\$\@.o] ;
  if ( make --version nothing 2> /dev/null | grep GNU > /dev/null ) ;  then
     foundGNUmake='yes, it is.' ;
     EXEC_DEPENDENCY='%: %.o' ;
  fi
  AC_MSG_RESULT($foundGNUmake)
  AC_SUBST(EXEC_DEPENDENCY)
])


################################################################
# Macro: ACX_MPI 
# Usage: ACX_MPI([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
# Authors: Steven G. Johnson and Julian C. Cummings
# Version: 2006-10-22
# Unmodified from the original; can be replaced with new version.
#
#      xref http://autoconf-archive.cryp.to/acx_mpi.html
#      Sets MPICC, MPILIBS output variable. 
#      If ACTION-IF-FOUND is not specified, default action defines HAVE_MPI.
#
AC_DEFUN([ACX_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

AC_LANG_CASE([C], [
        AC_REQUIRE([AC_PROG_CC])
        AC_ARG_VAR(MPICC,[MPI C compiler command])
        AC_CHECK_PROGS(MPICC, mpicc hcc mpxlc_r mpxlc mpcc cmpicc, $CC)
        acx_mpi_save_CC="$CC"
        CC="$MPICC"
        AC_SUBST(MPICC)
],
[C++], [
        AC_REQUIRE([AC_PROG_CXX])
        AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
        AC_CHECK_PROGS(MPICXX, mpic++ mpicxx mpiCC hcp mpxlC_r mpxlC mpCC cmpic++, $CXX)
        acx_mpi_save_CXX="$CXX"
        CXX="$MPICXX"
        AC_SUBST(MPICXX)
],
[Fortran 77], [
        AC_REQUIRE([AC_PROG_F77])
        AC_ARG_VAR(MPIF77,[MPI Fortran 77 compiler command])
        AC_CHECK_PROGS(MPIF77, mpif77 hf77 mpxlf_r mpxlf mpf77 cmpifc, $F77)
        acx_mpi_save_F77="$F77"
        F77="$MPIF77"
        AC_SUBST(MPIF77)
],
[Fortran], [
        AC_REQUIRE([AC_PROG_FC])
        AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
        AC_CHECK_PROGS(MPIFC, mpif90 mpxlf95_r mpxlf90_r mpxlf95 mpxlf90 mpf90 cmpif90c, $FC)
        acx_mpi_save_FC="$FC"
        FC="$MPIFC"
        AC_SUBST(MPIFC)
])

if test x = x"$MPILIBS"; then
        AC_LANG_CASE([C], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
                [C++], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
                [Fortran 77], [AC_MSG_CHECKING([for MPI_Init])
                        AC_LINK_IFELSE([AC_LANG_PROGRAM([],[      call MPI_Init])],[MPILIBS=" "
                                AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])],
                [Fortran], [AC_MSG_CHECKING([for MPI_Init])
                        AC_LINK_IFELSE([AC_LANG_PROGRAM([],[      call MPI_Init])],[MPILIBS=" "
                                AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])])
fi
AC_LANG_CASE([Fortran 77], [
        if test x = x"$MPILIBS"; then
                AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
        fi
        if test x = x"$MPILIBS"; then
                AC_CHECK_LIB(fmpich, MPI_Init, [MPILIBS="-lfmpich"])
        fi
],
[Fortran], [
        if test x = x"$MPILIBS"; then
                AC_CHECK_LIB(fmpi, MPI_Init, [MPILIBS="-lfmpi"])
        fi
        if test x = x"$MPILIBS"; then
                AC_CHECK_LIB(mpichf90, MPI_Init, [MPILIBS="-lmpichf90"])
        fi
])
if test x = x"$MPILIBS"; then
        AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
fi
if test x = x"$MPILIBS"; then
        AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
fi

dnl We have to use AC_TRY_COMPILE and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
AC_LANG_CASE([C], [if test x != x"$MPILIBS"; then
        AC_MSG_CHECKING([for mpi.h])
        AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
                AC_MSG_RESULT(no)])
fi],
[C++], [if test x != x"$MPILIBS"; then
        AC_MSG_CHECKING([for mpi.h])
        AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
                AC_MSG_RESULT(no)])
fi],
[Fortran 77], [if test x != x"$MPILIBS"; then
        AC_MSG_CHECKING([for mpif.h])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[      include 'mpif.h'])],[AC_MSG_RESULT(yes)], [MPILIBS=""
                AC_MSG_RESULT(no)])
fi],
[Fortran], [if test x != x"$MPILIBS"; then
        AC_MSG_CHECKING([for mpif.h])
        AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[      include 'mpif.h'])],[AC_MSG_RESULT(yes)], [MPILIBS=""
                AC_MSG_RESULT(no)])
fi])

AC_LANG_CASE([C], [CC="$acx_mpi_save_CC"],
        [C++], [CXX="$acx_mpi_save_CXX"],
        [Fortran 77], [F77="$acx_mpi_save_F77"],
        [Fortran], [FC="$acx_mpi_save_FC"])

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])dnl ACX_MPI



#################################################################
# Macro: ACX_PTHREAD
# Usage: ACX_PTHREAD([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]) 
# Authors:  Steven G. Johnson <stevenj@alum.mit.edu>
#           Alejandro Forero Cuervo <bachue@bachue.com>
# Version:  1.9 (2004/02/23)
# Source:   http://www.gnu.org/software/ac-archive/htmldoc/acx_pthread.html
# Everything below is verbatim from the archive. DO NOT MODIFY IT.
#
dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_pthread.html
dnl
AC_DEFUN([ACX_PTHREAD], [
AC_REQUIRE([AC_CANONICAL_HOST])
AC_LANG_SAVE
AC_LANG_C
acx_pthread_ok=no

# We used to check for pthread.h first, but this fails if pthread.h
# requires special compiler flags (e.g. on True64 or Sequent).
# It gets checked for in the link test anyway.

# First of all, check if the user has set any of the PTHREAD_LIBS,
# etcetera environment variables, and if threads linking works using
# them:
if test x"$PTHREAD_LIBS$PTHREAD_CFLAGS" != x; then
        save_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
        save_LIBS="$LIBS"
        LIBS="$PTHREAD_LIBS $LIBS"
        AC_MSG_CHECKING([for pthread_join in LIBS=$PTHREAD_LIBS with CFLAGS=$PTHREAD_CFLAGS])
        AC_TRY_LINK_FUNC(pthread_join, acx_pthread_ok=yes)
        AC_MSG_RESULT($acx_pthread_ok)
        if test x"$acx_pthread_ok" = xno; then
                PTHREAD_LIBS=""
                PTHREAD_CFLAGS=""
        fi
        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"
fi

# We must check for the threads library under a number of different
# names; the ordering is very important because some systems
# (e.g. DEC) have both -lpthread and -lpthreads, where one of the
# libraries is broken (non-POSIX).

# Create a list of thread flags to try.  Items starting with a "-" are
# C compiler flags, and other items are library names, except for "none"
# which indicates that we try without any flags at all, and "pthread-config"
# which is a program returning the flags for the Pth emulation library.

acx_pthread_flags="pthreads none -Kthread -kthread lthread -pthread -pthreads -mthreads pthread --thread-safe -mt pthread-config"

# The ordering *is* (sometimes) important.  Some notes on the
# individual items follow:

# pthreads: AIX (must check this before -lpthread)
# none: in case threads are in libc; should be tried before -Kthread and
#       other compiler flags to prevent continual compiler warnings
# -Kthread: Sequent (threads in libc, but -Kthread needed for pthread.h)
# -kthread: FreeBSD kernel threads (preferred to -pthread since SMP-able)
# lthread: LinuxThreads port on FreeBSD (also preferred to -pthread)
# -pthread: Linux/gcc (kernel threads), BSD/gcc (userland threads)
# -pthreads: Solaris/gcc
# -mthreads: Mingw32/gcc, Lynx/gcc
# -mt: Sun Workshop C (may only link SunOS threads [-lthread], but it
#      doesn't hurt to check since this sometimes defines pthreads too;
#      also defines -D_REENTRANT)
# pthread: Linux, etcetera
# --thread-safe: KAI C++
# pthread-config: use pthread-config program (for GNU Pth library)

case "${host_cpu}-${host_os}" in
        *solaris*)

        # On Solaris (at least, for some versions), libc contains stubbed
        # (non-functional) versions of the pthreads routines, so link-based
        # tests will erroneously succeed.  (We need to link with -pthread or
        # -lpthread.)  (The stubs are missing pthread_cleanup_push, or rather
        # a function called by this macro, so we could check for that, but
        # who knows whether they'll stub that too in a future libc.)  So,
        # we'll just look for -pthreads and -lpthread first:

        acx_pthread_flags="-pthread -pthreads pthread -mt $acx_pthread_flags"
        ;;
esac

if test x"$acx_pthread_ok" = xno; then
for flag in $acx_pthread_flags; do

        case $flag in
                none)
                AC_MSG_CHECKING([whether pthreads work without any flags])
                ;;

                -*)
                AC_MSG_CHECKING([whether pthreads work with $flag])
                PTHREAD_CFLAGS="$flag"
                ;;

		pthread-config)
		AC_CHECK_PROG(acx_pthread_config, pthread-config, yes, no)
		if test x"$acx_pthread_config" = xno; then continue; fi
		PTHREAD_CFLAGS="`pthread-config --cflags`"
		PTHREAD_LIBS="`pthread-config --ldflags` `pthread-config --libs`"
		;;

                *)
                AC_MSG_CHECKING([for the pthreads library -l$flag])
                PTHREAD_LIBS="-l$flag"
                ;;
        esac

        save_LIBS="$LIBS"
        save_CFLAGS="$CFLAGS"
        LIBS="$PTHREAD_LIBS $LIBS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"

        # Check for various functions.  We must include pthread.h,
        # since some functions may be macros.  (On the Sequent, we
        # need a special flag -Kthread to make this header compile.)
        # We check for pthread_join because it is in -lpthread on IRIX
        # while pthread_create is in libc.  We check for pthread_attr_init
        # due to DEC craziness with -lpthreads.  We check for
        # pthread_cleanup_push because it is one of the few pthread
        # functions on Solaris that doesn't have a non-functional libc stub.
        # We try pthread_create on general principles.
        AC_TRY_LINK([#include <pthread.h>],
                    [pthread_t th; pthread_join(th, 0);
                     pthread_attr_init(0); pthread_cleanup_push(0, 0);
                     pthread_create(0,0,0,0); pthread_cleanup_pop(0); ],
                    [acx_pthread_ok=yes])

        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"

        AC_MSG_RESULT($acx_pthread_ok)
        if test "x$acx_pthread_ok" = xyes; then
                break;
        fi

        PTHREAD_LIBS=""
        PTHREAD_CFLAGS=""
done
fi

# Various other checks:
if test "x$acx_pthread_ok" = xyes; then
        save_LIBS="$LIBS"
        LIBS="$PTHREAD_LIBS $LIBS"
        save_CFLAGS="$CFLAGS"
        CFLAGS="$CFLAGS $PTHREAD_CFLAGS"

        # Detect AIX lossage: threads are created detached by default
        # and the JOINABLE attribute has a nonstandard name (UNDETACHED).
        AC_MSG_CHECKING([for joinable pthread attribute])
        AC_TRY_LINK([#include <pthread.h>],
                    [int attr=PTHREAD_CREATE_JOINABLE;],
                    ok=PTHREAD_CREATE_JOINABLE, ok=unknown)
        if test x"$ok" = xunknown; then
                AC_TRY_LINK([#include <pthread.h>],
                            [int attr=PTHREAD_CREATE_UNDETACHED;],
                            ok=PTHREAD_CREATE_UNDETACHED, ok=unknown)
        fi
        if test x"$ok" != xPTHREAD_CREATE_JOINABLE; then
                AC_DEFINE(PTHREAD_CREATE_JOINABLE, $ok,
                          [Define to the necessary symbol if this constant
                           uses a non-standard name on your system.])
        fi
        AC_MSG_RESULT(${ok})
        if test x"$ok" = xunknown; then
                AC_MSG_WARN([we do not know how to create joinable pthreads])
        fi

        AC_MSG_CHECKING([if more special flags are required for pthreads])
        flag=no
        case "${host_cpu}-${host_os}" in
                *-aix* | *-freebsd*)     flag="-D_THREAD_SAFE";;
                *solaris* | *-osf* | *-hpux*) flag="-D_REENTRANT";;
        esac
        AC_MSG_RESULT(${flag})
        if test "x$flag" != xno; then
                PTHREAD_CFLAGS="$flag $PTHREAD_CFLAGS"
        fi

        LIBS="$save_LIBS"
        CFLAGS="$save_CFLAGS"

        # More AIX lossage: must compile with cc_r
        AC_CHECK_PROG(PTHREAD_CC, cc_r, cc_r, ${CC})
else
        PTHREAD_CC="$CC"
fi

AC_SUBST(PTHREAD_LIBS)
AC_SUBST(PTHREAD_CFLAGS)
AC_SUBST(PTHREAD_CC)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_pthread_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_PTHREAD,1,[Define if you have POSIX threads libraries and header files.]),[$1])
        :
else
        acx_pthread_ok=no
        $2
fi
AC_LANG_RESTORE
])dnl ACX_PTHREAD
#
# ACX_PTHREAD macro end.
# ****************************************************************
# ****************************************************************



#################################################################
# Macro: AX_CHECK_COMPILER_FLAGS
# Usage: AX_CHECK_COMPILER_FLAGS(FLAGS, [ACTION-SUCCESS], [ACTION-FAILURE])
# Authors:  Copyright (C) 2007 Steven G. Johnson <stevenj@alum.mit.edu>
#           Copyright (C) 2007 Matteo Frigo.
# Version:  2007-07-29
# Source:   http://autoconf-archive.cryp.to/ax_check_compiler_flags.html
#
# Check whether the given compiler FLAGS work with the current language's compiler, 
# or whether they give an error. (Warnings, however, are ignored.)
# ACTION-SUCCESS/ACTION-FAILURE are shell commands to execute on success/failure.
# 
# Everything below is verbatim from the archive. DO NOT MODIFY IT.
#
AC_DEFUN([AX_CHECK_COMPILER_FLAGS],
[AC_PREREQ(2.59) dnl for _AC_LANG_PREFIX
AC_MSG_CHECKING([whether _AC_LANG compiler accepts $1])
dnl Some hackery here since AC_CACHE_VAL can't handle a non-literal varname:
AS_LITERAL_IF([$1],
  [AC_CACHE_VAL(AS_TR_SH(ax_cv_[]_AC_LANG_ABBREV[]_flags_$1), [
      ax_save_FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
      _AC_LANG_PREFIX[]FLAGS="$1"
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
        AS_TR_SH(ax_cv_[]_AC_LANG_ABBREV[]_flags_$1)=yes,
        AS_TR_SH(ax_cv_[]_AC_LANG_ABBREV[]_flags_$1)=no)
      _AC_LANG_PREFIX[]FLAGS=$ax_save_FLAGS])],
  [ax_save_FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
   _AC_LANG_PREFIX[]FLAGS="$1"
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM()],
     eval AS_TR_SH(ax_cv_[]_AC_LANG_ABBREV[]_flags_$1)=yes,
     eval AS_TR_SH(ax_cv_[]_AC_LANG_ABBREV[]_flags_$1)=no)
   _AC_LANG_PREFIX[]FLAGS=$ax_save_FLAGS])
eval ax_check_compiler_flags=$AS_TR_SH(ax_cv_[]_AC_LANG_ABBREV[]_flags_$1)
AC_MSG_RESULT($ax_check_compiler_flags)
if test "x$ax_check_compiler_flags" = xyes; then
        m4_default([$2], :)
else
        m4_default([$3], :)
fi
])dnl AX_CHECK_COMPILER_FLAGS
#
# AX_CHECK_COMPILER_FLAGS macro end.
# ****************************************************************
# ****************************************************************




#################################################################
# Macro: AX_GCC_X86_CPUID
# Usage: AX_GCC_X86_CPUID(OP)
# Authors:  Copyright (C) 2007 Steven G. Johnson <stevenj@alum.mit.edu>
#           Copyright (C) 2007 Matteo Frigo
# Version:  2007-07-29
# Source:   http://autoconf-archive.cryp.to/ax_gcc_x86_cpuid.html
#
# Runs 'cpuid' with opcode 'OP'. 
# Sets cache variable ax_cv_gcc_x86_cpuid_OP to "eax:ebx:ecx:edx"
# where these are the registers set by 'cpuid'.
# If cpuid fails, variable is set to the string "unknown".
# This macro is required by AX_EXT; see below.
#
# Everything below is verbatim from the archive. DO NOT MODIFY IT.
#
AC_DEFUN([AX_GCC_X86_CPUID],
[AC_REQUIRE([AC_PROG_CC])
AC_LANG_PUSH([C])
AC_CACHE_CHECK(for x86 cpuid $1 output, ax_cv_gcc_x86_cpuid_$1,
 [AC_RUN_IFELSE([AC_LANG_PROGRAM([#include <stdio.h>], [
     int op = $1, eax, ebx, ecx, edx;
     FILE *f;
      __asm__("cpuid"
        : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx)
        : "a" (op));
     f = fopen("conftest_cpuid", "w"); if (!f) return 1;
     fprintf(f, "%x:%x:%x:%x\n", eax, ebx, ecx, edx);
     fclose(f);
     return 0;
])],
     [ax_cv_gcc_x86_cpuid_$1=`cat conftest_cpuid`; rm -f conftest_cpuid],
     [ax_cv_gcc_x86_cpuid_$1=unknown; rm -f conftest_cpuid],
     [ax_cv_gcc_x86_cpuid_$1=unknown])])
AC_LANG_POP([C])
])



#################################################################
# Macro: AX_EXT
# Usage: AX_EXT
# Author:   Copyright (C) 2007 Christophe Tournayre <turn3r@users.sourceforge.net>
# Version:  2007-11-28
# Source:   http://autoconf-archive.cryp.to/ax_ext.html
# Calls:    AC_SUBST(SIMD_CFLAGS)
# Defines:  HAVE_MMX / HAVE_SSE / HAVE_SSE2 / HAVE_SSE3 / HAVE_SSSE3          
#
# I commented out the bits that would define 
# HAVE_MMX / HAVE_SSE / HAVE_SSE3 / HAVE_SSSE3          
# and would have added -mmmx, etc to the SIMD_CFLAGS. 
# We only need -msse2
#
AC_DEFUN([AX_EXT],
[
  AC_REQUIRE([AX_GCC_X86_CPUID])

  AX_GCC_X86_CPUID(0x00000001)
  ecx=`echo $ax_cv_gcc_x86_cpuid_0x00000001 | cut -d ":" -f 3`
  edx=`echo $ax_cv_gcc_x86_cpuid_0x00000001 | cut -d ":" -f 4`

 AC_CACHE_CHECK([whether mmx is supported], [ax_have_mmx_ext],
  [
    ax_have_mmx_ext=no
    if test "$((0x$edx>>23&0x01))" = 1; then
      ax_have_mmx_ext=yes
    fi
  ])

 AC_CACHE_CHECK([whether sse is supported], [ax_have_sse_ext],
  [
    ax_have_sse_ext=no
    if test "$((0x$edx>>25&0x01))" = 1; then
      ax_have_sse_ext=yes
    fi
  ])

 AC_CACHE_CHECK([whether sse2 is supported], [ax_have_sse2_ext],
  [
    ax_have_sse2_ext=no
    if test "$((0x$edx>>26&0x01))" = 1; then
      ax_have_sse2_ext=yes
    fi
  ])

 AC_CACHE_CHECK([whether sse3 is supported], [ax_have_sse3_ext],
  [
    ax_have_sse3_ext=no
    if test "$((0x$ecx&0x01))" = 1; then
      ax_have_sse3_ext=yes
    fi
  ])

 AC_CACHE_CHECK([whether ssse3 is supported], [ax_have_ssse3_ext],
  [
    ax_have_ssse3_ext=no
    if test "$((0x$ecx>>9&0x01))" = 1; then
      ax_have_ssse3_ext=yes
    fi
  ])

#   if test "$ax_have_mmx_ext" = yes; then
#     AC_DEFINE(HAVE_MMX,,[Support mmx instructions])
#     AX_CHECK_COMPILER_FLAGS(-mmmx, SIMD_CFLAGS="$SIMD_CFLAGS -mmmx", [])
#   fi

#   if test "$ax_have_sse_ext" = yes; then
#     AC_DEFINE(HAVE_SSE,,[Support SSE (Streaming SIMD Extensions) instructions])
#     AX_CHECK_COMPILER_FLAGS(-msse, SIMD_CFLAGS="$SIMD_CFLAGS -msse", [])
#   fi

  if test "$ax_have_sse2_ext" = yes; then
    AC_DEFINE(HAVE_SSE2,,[Support SSE2 (Streaming SIMD Extensions 2) instructions])
    AX_CHECK_COMPILER_FLAGS(-msse2, SIMD_CFLAGS="$SIMD_CFLAGS -msse2", [])
  fi

#   if test "$ax_have_sse3_ext" = yes; then
#     AC_DEFINE(HAVE_SSE3,,[Support SSE3 (Streaming SIMD Extensions 3) instructions])
#     AX_CHECK_COMPILER_FLAGS(-msse3, SIMD_CFLAGS="$SIMD_CFLAGS -msse3", [])
#   fi

#   if test "$ax_have_ssse3_ext" = yes; then
#     AC_DEFINE(HAVE_SSSE3,,[Support SSSE3 (Supplemental Streaming SIMD Extensions 3) instructions])
#   fi

  AC_SUBST(SIMD_CFLAGS)
])
#
# AX_EXT macro end.
# ****************************************************************
# ****************************************************************



#################################################################
# Macro: AX_COMPILER_VENDOR
# Usage: AX_COMPILER_VENDOR
# Authors:  Copyright (C) 2007 Steven G. Johnson <stevenj@alum.mit.edu>
#           Copyright (C) 2007 Matteo Frigo
# Version:  2007-08-01
# Source:   http://autoconf-archive.cryp.to/ax_compiler_vendor.html
#
# Sets $ax_cv_c_compiler_vendor to gnu, intel, ibm, sun, hp, borland,
# comeau, dec, cray, kai, lcc, metrowerks, sgi, microsoft, watcom, etc.
#
# Everything below is verbatim from the archive. DO NOT MODIFY IT.
AC_DEFUN([AX_COMPILER_VENDOR],
[
AC_CACHE_CHECK([for _AC_LANG compiler vendor], ax_cv_[]_AC_LANG_ABBREV[]_compiler_vendor,
 [ax_cv_[]_AC_LANG_ABBREV[]_compiler_vendor=unknown
  # note: don't check for gcc first since some other compilers define __GNUC__
  for ventest in intel:__ICC,__ECC,__INTEL_COMPILER ibm:__xlc__,__xlC__,__IBMC__,__IBMCPP__ pathscale:__PATHCC__,__PATHSCALE__ gnu:__GNUC__ sun:__SUNPRO_C,__SUNPRO_CC hp:__HP_cc,__HP_aCC dec:__DECC,__DECCXX,__DECC_VER,__DECCXX_VER borland:__BORLANDC__,__TURBOC__ comeau:__COMO__ cray:_CRAYC kai:__KCC lcc:__LCC__ metrowerks:__MWERKS__ sgi:__sgi,sgi microsoft:_MSC_VER watcom:__WATCOMC__ portland:__PGI; do
    vencpp="defined("`echo $ventest | cut -d: -f2 | sed 's/,/) || defined(/g'`")"
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(,[
#if !($vencpp)
      thisisanerror;
#endif
])], [ax_cv_]_AC_LANG_ABBREV[_compiler_vendor=`echo $ventest | cut -d: -f1`; break])
  done
 ])
])
#
# AX_COMPILER_VENDOR macro end.
# ****************************************************************
# ****************************************************************



#################################################################
# Macro: AX_GCC_ARCHFLAG
# Usage: AX_GCC_ARCHFLAG([PORTABLE?], [ACTION-SUCCESS], [ACTION-FAILURE])
# Authors:  Copyright (C) 2007 Steven G. Johnson <stevenj@alum.mit.edu>
#           Copyright (C) 2007 Matteo Frigo
# Version:  2007-07-29
# Source:   http://autoconf-archive.cryp.to/ax_gcc_archflag.html
#
# This macro tries to guess the "native" arch corresponding to the
# target architecture for use with gcc's -march=arch or -mtune=arch
# flags. If found, the cache variable $ax_cv_gcc_archflag is set to this
# flag and ACTION-SUCCESS is executed; otherwise $ax_cv_gcc_archflag is
# is set to "unknown" and ACTION-FAILURE is executed. The default
# ACTION-SUCCESS is to add $ax_cv_gcc_archflag to the end of $CFLAGS.
#
# PORTABLE? should be either [yes] (default) or [no]. In the former
# case, the flag is set to -mtune (or equivalent) so that the
# architecture is only used for tuning, but the instruction set used is
# still portable. In the latter case, the flag is set to -march (or
# equivalent) so that architecture-specific instructions are enabled.
#
# The user can specify --with-gcc-arch=<arch> in order to override the
# macro's choice of architecture, or --without-gcc-arch to disable this.
#
# When cross-compiling, or if $CC is not gcc, then ACTION-FAILURE is
# called unless the user specified --with-gcc-arch manually.
#
# Everything below is verbatim from the archive. DO NOT MODIFY IT.
#
AC_DEFUN([AX_GCC_ARCHFLAG],
[AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_CANONICAL_HOST])

AC_ARG_WITH(gcc-arch, [AC_HELP_STRING([--with-gcc-arch=<arch>], [use architecture <arch> for gcc -march/-mtune, instead of guessing])],
        ax_gcc_arch=$withval, ax_gcc_arch=yes)

AC_MSG_CHECKING([for gcc architecture flag])
AC_MSG_RESULT([])
AC_CACHE_VAL(ax_cv_gcc_archflag,
[
ax_cv_gcc_archflag="unknown"

if test "$GCC" = yes; then

if test "x$ax_gcc_arch" = xyes; then
ax_gcc_arch=""
if test "$cross_compiling" = no; then
case $host_cpu in
  i[[3456]]86*|x86_64*) # use cpuid codes, in part from x86info-1.7 by D. Jones
     AX_GCC_X86_CPUID(0)
     AX_GCC_X86_CPUID(1)
     case $ax_cv_gcc_x86_cpuid_0 in
       *:756e6547:*:*) # Intel
          case $ax_cv_gcc_x86_cpuid_1 in
            *5[[48]]?:*:*:*) ax_gcc_arch="pentium-mmx pentium" ;;
            *5??:*:*:*) ax_gcc_arch=pentium ;;
            *6[[3456]]?:*:*:*) ax_gcc_arch="pentium2 pentiumpro" ;;
            *6a?:*[[01]]:*:*) ax_gcc_arch="pentium2 pentiumpro" ;;
            *6a?:*[[234]]:*:*) ax_gcc_arch="pentium3 pentiumpro" ;;
            *6[[9d]]?:*:*:*) ax_gcc_arch="pentium-m pentium3 pentiumpro" ;;
            *6[[78b]]?:*:*:*) ax_gcc_arch="pentium3 pentiumpro" ;;
            *6??:*:*:*) ax_gcc_arch=pentiumpro ;;
            *f3[[347]]:*:*:*|*f4[1347]:*:*:*)
                case $host_cpu in
                  x86_64*) ax_gcc_arch="nocona pentium4 pentiumpro" ;;
                  *) ax_gcc_arch="prescott pentium4 pentiumpro" ;;
                esac ;;
            *f??:*:*:*) ax_gcc_arch="pentium4 pentiumpro";;
          esac ;;
       *:68747541:*:*) # AMD
          case $ax_cv_gcc_x86_cpuid_1 in
            *5[[67]]?:*:*:*) ax_gcc_arch=k6 ;;
            *5[[8d]]?:*:*:*) ax_gcc_arch="k6-2 k6" ;;
            *5[[9]]?:*:*:*) ax_gcc_arch="k6-3 k6" ;;
            *60?:*:*:*) ax_gcc_arch=k7 ;;
            *6[[12]]?:*:*:*) ax_gcc_arch="athlon k7" ;;
            *6[[34]]?:*:*:*) ax_gcc_arch="athlon-tbird k7" ;;
            *67?:*:*:*) ax_gcc_arch="athlon-4 athlon k7" ;;
            *6[[68a]]?:*:*:*)
               AX_GCC_X86_CPUID(0x80000006) # L2 cache size
               case $ax_cv_gcc_x86_cpuid_0x80000006 in
                 *:*:*[[1-9a-f]]??????:*) # (L2 = ecx >> 16) >= 256
                        ax_gcc_arch="athlon-xp athlon-4 athlon k7" ;;
                 *) ax_gcc_arch="athlon-4 athlon k7" ;;
               esac ;;
            *f[[4cef8b]]?:*:*:*) ax_gcc_arch="athlon64 k8" ;;
            *f5?:*:*:*) ax_gcc_arch="opteron k8" ;;
            *f7?:*:*:*) ax_gcc_arch="athlon-fx opteron k8" ;;
            *f??:*:*:*) ax_gcc_arch="k8" ;;
          esac ;;
        *:746e6543:*:*) # IDT
           case $ax_cv_gcc_x86_cpuid_1 in
             *54?:*:*:*) ax_gcc_arch=winchip-c6 ;;
             *58?:*:*:*) ax_gcc_arch=winchip2 ;;
             *6[[78]]?:*:*:*) ax_gcc_arch=c3 ;;
             *69?:*:*:*) ax_gcc_arch="c3-2 c3" ;;
           esac ;;
     esac
     if test x"$ax_gcc_arch" = x; then # fallback
        case $host_cpu in
          i586*) ax_gcc_arch=pentium ;;
          i686*) ax_gcc_arch=pentiumpro ;;
        esac
     fi
     ;;

  sparc*)
     AC_PATH_PROG([PRTDIAG], [prtdiag], [prtdiag], [$PATH:/usr/platform/`uname -i`/sbin/:/usr/platform/`uname -m`/sbin/])
     cputype=`(((grep cpu /proc/cpuinfo | cut -d: -f2) ; ($PRTDIAG -v |grep -i sparc) ; grep -i cpu /var/run/dmesg.boot ) | head -n 1) 2> /dev/null`
     cputype=`echo "$cputype" | tr -d ' -' |tr $as_cr_LETTERS $as_cr_letters`
     case $cputype in
         *ultrasparciv*) ax_gcc_arch="ultrasparc4 ultrasparc3 ultrasparc v9" ;;
         *ultrasparciii*) ax_gcc_arch="ultrasparc3 ultrasparc v9" ;;
         *ultrasparc*) ax_gcc_arch="ultrasparc v9" ;;
         *supersparc*|*tms390z5[[05]]*) ax_gcc_arch="supersparc v8" ;;
         *hypersparc*|*rt62[[056]]*) ax_gcc_arch="hypersparc v8" ;;
         *cypress*) ax_gcc_arch=cypress ;;
     esac ;;

  alphaev5) ax_gcc_arch=ev5 ;;
  alphaev56) ax_gcc_arch=ev56 ;;
  alphapca56) ax_gcc_arch="pca56 ev56" ;;
  alphapca57) ax_gcc_arch="pca57 pca56 ev56" ;;
  alphaev6) ax_gcc_arch=ev6 ;;
  alphaev67) ax_gcc_arch=ev67 ;;
  alphaev68) ax_gcc_arch="ev68 ev67" ;;
  alphaev69) ax_gcc_arch="ev69 ev68 ev67" ;;
  alphaev7) ax_gcc_arch="ev7 ev69 ev68 ev67" ;;
  alphaev79) ax_gcc_arch="ev79 ev7 ev69 ev68 ev67" ;;

  powerpc*)
     cputype=`((grep cpu /proc/cpuinfo | head -n 1 | cut -d: -f2 | cut -d, -f1 | sed 's/ //g') ; /usr/bin/machine ; /bin/machine; grep CPU /var/run/dmesg.boot | head -n 1 | cut -d" " -f2) 2> /dev/null`
     cputype=`echo $cputype | sed -e 's/ppc//g;s/ *//g'`
     case $cputype in
       *750*) ax_gcc_arch="750 G3" ;;
       *740[[0-9]]*) ax_gcc_arch="$cputype 7400 G4" ;;
       *74[[4-5]][[0-9]]*) ax_gcc_arch="$cputype 7450 G4" ;;
       *74[[0-9]][[0-9]]*) ax_gcc_arch="$cputype G4" ;;
       *970*) ax_gcc_arch="970 G5 power4";;
       *POWER4*|*power4*|*gq*) ax_gcc_arch="power4 970";;
       *POWER5*|*power5*|*gr*|*gs*) ax_gcc_arch="power5 power4 970";;
       603ev|8240) ax_gcc_arch="$cputype 603e 603";;
       *) ax_gcc_arch=$cputype ;;
     esac
     ax_gcc_arch="$ax_gcc_arch powerpc"
     ;;
esac
fi # not cross-compiling
fi # guess arch

if test "x$ax_gcc_arch" != x -a "x$ax_gcc_arch" != xno; then
for arch in $ax_gcc_arch; do
  if test "x[]m4_default([$1],yes)" = xyes; then # if we require portable code
    flags="-mtune=$arch"
    # -mcpu=$arch and m$arch generate nonportable code on every arch except
    # x86.  And some other arches (e.g. Alpha) don't accept -mtune.  Grrr.
    case $host_cpu in i*86|x86_64*) flags="$flags -mcpu=$arch -m$arch";; esac
  else
    flags="-march=$arch -mcpu=$arch -m$arch"
  fi
  for flag in $flags; do
    AX_CHECK_COMPILER_FLAGS($flag, [ax_cv_gcc_archflag=$flag; break])
  done
  test "x$ax_cv_gcc_archflag" = xunknown || break
done
fi

fi # $GCC=yes
])
AC_MSG_CHECKING([for gcc architecture flag])
AC_MSG_RESULT($ax_cv_gcc_archflag)
if test "x$ax_cv_gcc_archflag" = xunknown; then
  m4_default([$3],:)
else
  m4_default([$2], [CFLAGS="$CFLAGS $ax_cv_gcc_archflag"])
fi
])
#
# AX_GCC_ARCHFLAG macro end.
# ****************************************************************
# ****************************************************************


#################################################################
# Macro: AX_CC_MAXOPT
# Usage: AX_CC_MAXOPT
# Authors:  Copyright (C) 2007 Steven G. Johnson <stevenj@alum.mit.edu>
#           Copyright (C) 2007 Matteo Frigo
# Version:  2007-07-29
# Source:   http://autoconf-archive.cryp.to/ax_cc_maxopt.html
#
# Try to turn on "good" C optimization flags for various compilers and
# architectures, for some definition of "good". 
#
# The user can override the flags by setting the CFLAGS environment
# variable.  The user can also specify --enable-portable-binary in
# order to disable any optimization flags that might result in a
# binary that only runs on the host architecture.
#
# Note also that the flags assume that ANSI C aliasing rules are
# followed by the code (e.g. for gcc's -fstrict-aliasing), and that
# floating-point computations can be re-ordered as needed.
#
# SRE: I've made modifications as follows.
#  - HMMER relies on IEEE754-compliant math. Don't enable
#    any options that break compliance; for example, gcc -ffast-math
#
AC_DEFUN([AX_CC_MAXOPT],
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AX_COMPILER_VENDOR])
AC_REQUIRE([AC_CANONICAL_HOST])

AC_ARG_ENABLE(portable-binary, [AC_HELP_STRING([--enable-portable-binary], [disable compiler optimizations that would produce unportable binaries])],
        acx_maxopt_portable=$withval, acx_maxopt_portable=no)

# Try to determine "good" native compiler flags if none specified via CFLAGS
if test "$ac_test_CFLAGS" != "set"; then
  CFLAGS=""
  case $ax_cv_c_compiler_vendor in
    dec) CFLAGS="-newc -w0 -O5 -ansi_alias -ansi_args -tune host"
#        CFLAGS="-newc -w0 -O5 -ansi_alias -ansi_args -fp_reorder -tune host"
         if test "x$acx_maxopt_portable" = xno; then
           CFLAGS="$CFLAGS -arch host"
         fi;;

    sun) CFLAGS="-native -xO5 -dalign"
#        CFLAGS="-native -fast -xO5 -dalign"
         if test "x$acx_maxopt_portable" = xyes; then
           CFLAGS="$CFLAGS -xarch=generic"
         fi;;

    hp)  CFLAGS="+Oall +Optrs_ansi +DSnative"
         if test "x$acx_maxopt_portable" = xyes; then
           CFLAGS="$CFLAGS +DAportable"
         fi;;

    ibm) if test "x$acx_maxopt_portable" = xno; then
           xlc_opt="-qarch=auto -qtune=auto"
         else
           xlc_opt="-qtune=auto"
         fi
         AX_CHECK_COMPILER_FLAGS($xlc_opt,
                CFLAGS="-O3 -qansialias -w $xlc_opt",
               [CFLAGS="-O3 -qansialias -w"
                echo "******************************************************"
                echo "*  You seem to have the IBM  C compiler.  It is      *"
                echo "*  recommended for best performance that you use:    *"
                echo "*                                                    *"
                echo "*    CFLAGS=-O3 -qarch=xxx -qtune=xxx -qansialias -w *"
                echo "*                      ^^^        ^^^                *"
                echo "*  where xxx is pwr2, pwr3, 604, or whatever kind of *"
                echo "*  CPU you have.  (Set the CFLAGS environment var.   *"
                echo "*  and re-run configure.)  For more info, man cc.    *"
                echo "******************************************************"])
         ;;

    intel) CFLAGS="-O3 -ansi_alias"
        if test "x$acx_maxopt_portable" = xno; then
          icc_archflag=unknown
          icc_flags=""
          case $host_cpu in
            i686*|x86_64*)
              # icc accepts gcc assembly syntax, so these should work:
              AX_GCC_X86_CPUID(0)
              AX_GCC_X86_CPUID(1)
              case $ax_cv_gcc_x86_cpuid_0 in # see AX_GCC_ARCHFLAG
                *:756e6547:*:*) # Intel
                  case $ax_cv_gcc_x86_cpuid_1 in
                    *6a?:*[[234]]:*:*|*6[[789b]]?:*:*:*) icc_flags="-xK";;
                    *f3[[347]]:*:*:*|*f4[1347]:*:*:*) icc_flags="-xP -xN -xW -xK";;
                    *f??:*:*:*) icc_flags="-xN -xW -xK";;
                  esac ;;
              esac ;;
          esac
          if test "x$icc_flags" != x; then
            for flag in $icc_flags; do
              AX_CHECK_COMPILER_FLAGS($flag, [icc_archflag=$flag; break])
            done
          fi
          AC_MSG_CHECKING([for icc architecture flag])
          AC_MSG_RESULT($icc_archflag)
          if test "x$icc_archflag" != xunknown; then
            CFLAGS="$CFLAGS $icc_archflag"
          fi
        fi
        ;;

    gnu)
     # default optimization flags for gcc on all systems
     CFLAGS="-O3 -fomit-frame-pointer"

     # -malign-double for x86 systems
     AX_CHECK_COMPILER_FLAGS(-malign-double, CFLAGS="$CFLAGS -malign-double")

     #  -fstrict-aliasing for gcc-2.95+
     AX_CHECK_COMPILER_FLAGS(-fstrict-aliasing,
        CFLAGS="$CFLAGS -fstrict-aliasing")

     # note that we enable "unsafe" fp optimization with other compilers, too
     # SRE: no, that's a bad idea, don't use this
#     AX_CHECK_COMPILER_FLAGS(-ffast-math, CFLAGS="$CFLAGS -ffast-math")

     AX_GCC_ARCHFLAG($acx_maxopt_portable)
     ;;
  esac

  if test -z "$CFLAGS"; then
        echo ""
        echo "********************************************************"
        echo "* WARNING: Don't know the best CFLAGS for this system  *"
        echo "* Use ./configure CFLAGS=... to specify your own flags *"
        echo "* (otherwise, a default of CFLAGS=-O3 will be used)    *"
        echo "********************************************************"
        echo ""
        CFLAGS="-O3"
  fi

  AX_CHECK_COMPILER_FLAGS($CFLAGS, [], [
        echo ""
        echo "********************************************************"
        echo "* WARNING: The guessed CFLAGS don't seem to work with  *"
        echo "* your compiler.                                       *"
        echo "* Use ./configure CFLAGS=... to specify your own flags *"
        echo "********************************************************"
        echo ""
        CFLAGS=""
  ])

fi
])
#
# AX_CC_MAXOPT macro end.
# ****************************************************************
# ****************************************************************
