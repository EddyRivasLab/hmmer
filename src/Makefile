SUBDIRS = base\
	  build\
          daemon\
          dp_reference\
	  dp_sparse\
	  dp_vector\
	  misc\
	  sandbox \
          search\
	  programs \
          utilities \
          hardware\
	  experiments

HDRS =  

OBJS = 

STATS = 

BENCHMARKS = 

UTESTS =

ITESTS = 

EXAMPLES = 



# VPATH and shell configuration
top_srcdir   = ..
top_builddir = ..
srcdir       = .

SHELL        = /bin/sh

# Installation targets
#
prefix      = /usr/local
exec_prefix = ${prefix}
datarootdir = ${prefix}/share
bindir      = ${exec_prefix}/bin
libdir      = ${exec_prefix}/lib
includedir  = ${prefix}/include
mandir      = ${datarootdir}/man
docdir      = ${datarootdir}/doc/${PACKAGE_TARNAME}
pdfdir      = ${docdir}
mandir      = ${datarootdir}/man
man1dir     = ${mandir}/man1
man1ext     = .1

CC        = mpicc
CFLAGS    = -O3 -g -pthread -fPIC
SIMDFLAGS = @SIMD_CFLAGS@
CPPFLAGS  = 
LDFLAGS   = 
DEFS      = -DHAVE_CONFIG_H
LIBS      = -lhmmer -leasel -ldivsufsort      -lm

AR        = /usr/bin/ar 
RANLIB    = ranlib
INSTALL   = /usr/bin/install -c

ESLDIR    = lib/easel
SADIR     = lib/libdivsufsort

MYLIBDIRS = -L${top_builddir}/${ESLDIR} -L${top_builddir}/${SADIR} -L.
MYINCDIRS = -I${top_builddir}/${ESLDIR} \
            -I${top_builddir}/${SADIR} \
	    -I. \
            -I${top_srcdir}/${ESLDIR} \
	    -I${top_srcdir}/${SADIR} \
	    -I${top_srcdir}


# beautification magic stolen from git
QUIET_SUBDIR0 = +${MAKE} -C #space separator after -c
QUIET_SUBDIR1 = 
ifndef V
	QUIET_CC      = @echo '    ' CC $@;
	QUIET_GEN     = @echo '    ' GEN $@;
	QUIET_AR      = @echo '    ' AR $@;
	QUIET_SUBDIR0 = +@subdir=
	QUIET_SUBDIR1 = ; echo '    ' SUBDIR $$subdir; \
		        ${MAKE} -C $$subdir
	QUIET         = @
endif

.PHONY: all dev tests check install uninstall distclean clean TAGS tags-append

all: libhmmer.a 
	${QUIET}+set -e; for subdir in ${SUBDIRS}; do \
	   echo '    ' SUBDIR $$subdir; \
	   ${MAKE} -C $$subdir all; \
	done

dev: libhmmer.a
	${QUIET}+set -e; for subdir in ${SUBDIRS}; do \
	   echo '    ' SUBDIR $$subdir; \
	   ${MAKE} -C $$subdir dev; \
	done

tests: libhmmer.a
	${QUIET}+set -e; for subdir in ${SUBDIRS}; do \
	   echo '    ' SUBDIR $$subdir; \
	   ${MAKE} -C $$subdir tests; \
	done

check: libhmmer.a
	${QUIET}+set -e; for subdir in ${SUBDIRS}; do \
	   echo '    ' SUBDIR $$subdir; \
	   ${MAKE} -C $$subdir check; \
	done

# We need to do one pass thru the subdirs to build libhmmer.a before
# we start trying to build any executables. That's why we make
# libhmmer.a a target, and force building the libhmmer-{subdir}.stamp
# target in each subdir.
.FORCE:
libhmmer.a: .FORCE
	${QUIET}+set -e; for subdir in ${SUBDIRS}; do \
	   echo '    ' SUBDIR $$subdir; \
	   ${MAKE} -C $$subdir libhmmer-$$subdir.stamp; \
	done


install: libhmmer.a
	${QUIET}+set -e; for subdir in ${SUBDIRS}; do \
	   echo '    ' SUBDIR $$subdir; \
	   ${MAKE} -C $$subdir install ;\
	done

uninstall:
	${QUIET}+set -e; for subdir in ${SUBDIRS}; do \
	   echo '    ' SUBDIR $$subdir ;\
	   ${MAKE} -C $$subdir uninstall ;\
	done

distclean: clean
	${QUIET}+set -e; for subdir in ${SUBDIRS}; do \
	   echo '    ' SUBDIR $$subdir; \
	   ${MAKE} -C $$subdir distclean; \
	done
	-rm -f Makefile Makefile-subdirs.mk p7_config.h

clean:
	${QUIET}+set -e; for subdir in ${SUBDIRS}; do \
	   echo '    ' SUBDIR $$subdir; \
	   ${MAKE} -C $$subdir clean; \
	done
	-rm -f *.o *~ Makefile.bak core ${PROGS} TAGS gmon.out
	-rm -f libhmmer.a 

tags-append:
	${QUIET}+set -e; for subdir in ${SUBDIRS}; do \
	   echo '    ' SUBDIR $$subdir; \
	   ${MAKE} -C $$subdir tags-append; \
	done
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/*.h ${srcdir}/*.in


################################################################
# @LICENSE@
#
# SVN $URL$
# SVN $Id$
################################################################
