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

