f# Process this file with autoconf to produce a configure script.
# SRE, Mon Oct  5 14:55:45 1998
# CVS $Id$
# 
# Hint: See autoconf macro archive at 
#     http://www.gnu.org/software/ac-archive/

# Initialize w/ AC_INIT(<package> <version> <bug-report> <tarname>)
# 
# The four AC_INIT args set the following output variables and preprocessor symbols:
#     PACKAGE_NAME      <package>     e.g. "HMMER"
#     PACKAGE_VERSION   <version>     e.g. "2.2h"
#     PACKAGE_BUGREPORT <bug-report>  e.g. "eddy@genetics.wustl.edu"
#     PACKAGE_TARNAME   <tarname>     e.g. "hmmer"
# From them, it derives one more:
#     PACKAGE_STRING    <package> <version>, e.g. "HMMER 2.2h"
# We also define additional variables:
#     PACKAGE_RELCODE     unique CVS tag without .'s: e.g. "hmmer2_2i"
#     PACKAGE_DATE        release date: e.g. "January 2003"
#     PACKAGE_COPYRIGHT   one-line copyright string
#     PACKAGE_LICENSE     one-line license string
#     PACKAGE_LICENSETAG  which license to bundle from Licenses/ subdirectory.
#
# The preprocessor symbols get picked up in config.h.
# The output variables are used in Makefiles.
#
AC_PREREQ(2.57)
AC_INIT(HMMER, 2.4dev1, eddy@genetics.wustl.edu, hmmer)
AC_MSG_NOTICE([configuring HMMER for your system.])

PACKAGE_RELCODE="hmmer2_4dev1"
PACKAGE_DATE="June 2004"
PACKAGE_COPYRIGHT="Copyright (C) 1992-2004 HHMI/Washington University School of Medicine"
PACKAGE_LICENSE="Freely distributed under the GNU General Public License (GPL)"
PACKAGE_LICENSETAG=gnu




# For a normal release, you don't need to edit below this line.
#
AC_SUBST(PACKAGE_RELCODE)
AC_SUBST(PACKAGE_DATE)
AC_SUBST(PACKAGE_COPYRIGHT)
AC_SUBST(PACKAGE_LICENSE)
AC_SUBST(PACKAGE_LICENSETAG)

AC_DEFINE_UNQUOTED(PACKAGE_DATE, "$PACKAGE_DATE")
AC_DEFINE_UNQUOTED(PACKAGE_COPYRIGHT, "$PACKAGE_COPYRIGHT")
AC_DEFINE_UNQUOTED(PACKAGE_LICENSE, "$PACKAGE_LICENSE")


# Check if we're a binary distribution. 
# The trigger for this is the existence of the binaries/ subdirectory,
# and the nonexistence of a src/ subdirectory.
# 
# If we are a binary dist, we only do a minimal configuration, setting up
# version information and installation path names.
# This makes ./configure independent of having a working
# C compiler.
#
# Watch out: this is implemented as a sneaky if statement that wraps the
# majority of the configure.ac M4 code.
#
if ((! test -d binaries) && test -d src); then 
AC_MSG_NOTICE([This is a full source distribution - using full config])

# Checks for programs.
# Sets output variables CC, CFLAGS, and RANLIB.
# AR might be hidden away on Solaris; provide some hints in the PATH.
#
AC_PROG_CC
if test "$CC" = xlc ; then	# IBM's xlc?
  CFLAGS=$CFLAGS" -O3 -qarch=auto -qtune=auto"
fi
AC_PROG_RANLIB
AC_PATH_PROG([AR], [ar], [:], [$PATH:/usr/ccs/bin:/usr/xpg4/bin])


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



# GNU make check; format of makefile dependency lines for executables.
# original from John Darrington <j.darrington@elvis.murdoch.edu.au> 
# w/ heavy modifications.
#
# We need this because GNU make and SYSV make use different systems
# specifying variables for dependencies: $$@ in sysv, %: %.o in GNU.
# Would love to hear a better way of doing this.
# 
# I use two different conventions in my Makefiles. Sometimes 
# executable "foo" has a file "foo.c" - this is the HMMER convention.
# Sometimes executable "foo" has a file "foo_main.c" - this is
# the SQUID convention. The configure script sets the
# EXEC_DEPENDENCY appropriately: here, HMMER style.
#
# This creates a function CHECK_GNU_MAKE, which we immediately call.
# It sets an output variable EXEC_DEPENDENCY. 
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
CHECK_GNU_MAKE

#################################################################################
# Now we're into our "optional features"
#################################################################################


# --enable-altivec         - enable Erik Lindahl's Altivec code for PowerPC
#
# Sets the ALTIVEC preprocessor variable; adds -faltivec to the CFLAGS.
#
AC_ARG_ENABLE(altivec,
[  --enable-altivec        enable Altivec-specific code for PowerPC],
[case $enable_altivec in
   yes)  AC_MSG_NOTICE([enabled Altivec optimization for PowerPC (Erik Lindahl)])
         AC_DEFINE(ALTIVEC)
	# We support altivec both on Apple platforms (Lindahl's contribs)
        # and on IBM PowerPC platforms (IBM/Yinhe Cheng's contribs).
	# Three cases: 
	#  1. we're using cc/gcc on an IBM PowerPC;
        #  2. we're using OS/X's cc on an Apple;
        #  3. we're using IBM's xlc on an IBM PowerPC.
	 if test "$CC" != xlc; then 
	   if test "$ac_cv_build" = powerpc64-unknown-linux-gnu; then 
	     altivecOP="-maltivec -mabi=altivec"  # IBM gcc
	   else
	     altivecOP="-faltivec -mcpu=G5 -mtune=G5" # Apple cc
	   fi
	 else
	   altivecOP="-qaltivec"  # IBM xlc
	 fi

         CFLAGS="$CFLAGS $altivecOP" 
         ;;
    no)  AC_MSG_NOTICE([Altivec/PowerPC optimization disabled])
         ;;
     *)  echo "Ignoring unknown argument to --enable-altivec: $enable_altivec"
	 ;;
esac])
#
# Now check whether the compiler supports altivec extensions, if enabled.
# In theory we could do this in the source, but this way we can give
# the user an explicit error message telling him/her to upgrade the compiler.
# (M4 code from Erik Lindahl, Stanford University; modified to include
# IBM PowerPC xlc tests by Yinhe Cheng, IBM).
#
if test x"$enable_altivec" = xyes; then
	if test "$ac_cv_build" = powerpc64-unknown-linux-gnu; then
	   cat >> confdefs.h << \_ACEOF
#include "altivec.h"
_ACEOF
        fi

        AC_MSG_CHECKING([whether the compiler supports altivec extensions])
############################################################################
#### IBM PPC xlc altivec compile test: (Yinhe Cheng, IBM)
############################################################################
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
############################################################################
#### else, Apple OS/X altivec compile test: (Erik Lindahl, Stanford)
############################################################################
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




# --enable-ccmalloc      - configure for compiling in ccmalloc memory checks
#
# The CC output variable is modified.
# Requires ccmalloc version >= 0.4.0 because of the --no-wrapper option.
#
AC_ARG_ENABLE(ccmalloc,
[  --enable-ccmalloc       turn ccmalloc memory debugging on ],
[case $enable_ccmalloc in
   yes)  AC_MSG_NOTICE([enabled ccmalloc memory debugging])
         CC="ccmalloc --no-wrapper $CC" 
         ;;
    no)  AC_MSG_NOTICE([ccmalloc memory debug/tracing disabled])
         ;;
     *)  echo "Ignoring unknown argument to --enable-ccmalloc: $enable_ccmalloc"
	 ;;
esac])


# --enable-debugging=x    - set debugging level to <x> (1-3)
# 
# At all levels, including 0, replaces CFLAGS w/ "-g -Wall" (so it assumes gcc).
# Sets the DEBUGLEVEL preprocessor symbol to <x>
#
AC_ARG_ENABLE(debugging,
[  --enable-debugging      set CFLAGS for debugging
  --enable-debugging=x    also set diagnostics level to <x> (1-3) 
                          (where 3 = most verbose)],
[ case $enable_debugging in
   yes)  AC_MSG_NOTICE([enabled debugging diagnostics level 0 (CFLAGS only, no verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(DEBUGLEVEL, 0)
         ;;
     1)  AC_MSG_NOTICE([enabled debugging diagnostics level 1 (low verbosity)])
         CFLAGS="-g -Wall"
         AC_DEFINE(DEBUGLEVEL, 1)
         ;;
     2)  AC_MSG_NOTICE([enabled debugging diagnostics level 2 (moderate verbosity)])
         CFLAGS="-g -Wall"
	 AC_DEFINE(DEBUGLEVEL, 2)
         ;;
     3)  AC_MSG_NOTICE([enabled debugging diagnostics level 3 (high verbosity)])
         CFLAGS="-g -Wall"
	 AC_DEFINE(DEBUGLEVEL, 3)
         ;;
    no)  AC_MSG_NOTICE([debugging diagnostics disabled])
         ;;
     *)  echo "Ignoring unknown argument to --enable-debugging: $enable_debugging"
	 ;;
esac])


# --enable-lfs           Large File Summit (LFS) support for >2GB files
# See: http://ftp.sas.com/standards/large.file/x_open.20Mar96.html
#
AC_ARG_ENABLE(lfs, 
[  --enable-lfs            enable LFS, Large File Support],
[case $enable_lfs in
   yes) AC_MSG_NOTICE([configured for optional LFS, large file support])
	AC_DEFINE(_LARGEFILE_SOURCE)
	AC_DEFINE(_LARGEFILE64_SOURCE)
	AC_DEFINE(_FILE_OFFSET_BITS, 64)
        ;;
   no)  ;;
    *)  echo "Ignoring unknown argument to --enable-lfs: $enable_lfs"
        ;;
esac])



# --enable-pvm       Enable Parallel Virtual Machine (PVM) support
# 
# Sets PVMLIBDIR, PVMINCDIR, PCMPROGS, PVMLIBS output variables
# Sets HMMER_PVM preprocessor variable.
#
AC_ARG_ENABLE(pvm, 
[  --enable-pvm            enable PVM, Parallel Virtual Machine],
[case $enable_pvm in
   yes)  AC_MSG_NOTICE([enabled optional PVM (Parallel Virtual Machine) support])
 	 PVMLIBDIR="-L${PVM_ROOT}/lib/${PVM_ARCH}"
	 PVMINCDIR="-I${PVM_ROOT}/include"
	 PVMPROGS="hmmcalibrate-pvm hmmpfam-pvm hmmsearch-pvm"
	 PVMLIBS="-lpvm3"
	 AC_DEFINE(HMMER_PVM)
	 ;;
   no)   AC_MSG_NOTICE([PVM (Parallel Virtual machine) support disabled])
         ;;
   *)    echo "Ignoring unknown argument to --enable-pvm: $enable_pvm"
	 ;;
esac])
AC_SUBST(PVMLIBDIR)
AC_SUBST(PVMINCDIR)
AC_SUBST(PVMPROGS)
AC_SUBST(PVMLIBS)


# --enable-threads      Enable POSIX multithreading
#
# Uses ACX_PTHREAD macro from the GNU macro archive.
# ****************************************************************
# ****************************************************************
# 
# The ACX_PTHREAD macro.
# Synopsis: ACX_PTHREAD([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]) 
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
# ACX_PTHREAD macro end.
# ****************************************************************
# ****************************************************************

#
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


# Check for incompatible optional features
#
if test x"$enable_threads" = xyes; then
  if test x"$enable_ccmalloc" = xyes; then
    AC_MSG_ERROR([ccmalloc probably not threadsafe; do not enable both ccmalloc and threads])
  fi
fi


# Although we configure our fixed-size integers in squid,
# we do need some of that configuration code here, because
# we might (on Solaris, for instance) need to include
# -lsocket to the linker to get ntohs() and friends. These
# define preprocessor variables HAVE_NTOHS, HAVE_NTOHL,
# HAVE_HTONS, and HAVE_HTONL.
#
AC_CHECK_FUNCS(ntohs, , AC_CHECK_LIB(socket, ntohs))
AC_CHECK_FUNCS(ntohl, , AC_CHECK_LIB(socket, ntohl))
AC_CHECK_FUNCS(htons, , AC_CHECK_LIB(socket, htons))
AC_CHECK_FUNCS(htonl, , AC_CHECK_LIB(socket, htonl))

# Configure the squid and easel libraries
#
AC_CONFIG_SUBDIRS(squid)
AC_CONFIG_SUBDIRS(easel)

# Write the config.h header
#
AC_CONFIG_HEADERS(src/config.h)

# Write out Makefiles
#
AC_CONFIG_FILES(src/Makefile)
AC_CONFIG_FILES(testsuite/Makefile)

# remember that sneaky if statement that wraps the source
# configuration? well, it ends here.
else
  AC_MSG_NOTICE([This is a precompiled binary distribution - using abbreviated config])
fi

AC_CONFIG_FILES(Makefile)
AC_OUTPUT



