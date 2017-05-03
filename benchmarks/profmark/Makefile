PROGS        = create-profmark\
	       rocplot

SHELL        = /bin/sh
top_srcdir   = ../..
top_builddir = ../..
srcdir       = .


prefix       = /usr/local
exec_prefix  = ${prefix}
datarootdir  = ${prefix}/share
bindir       = ${exec_prefix}/bin
libdir       = ${exec_prefix}/lib
includedir   = ${prefix}/include

CC           = mpicc
CFLAGS       = -g -Wall -pthread -fPIC
#SIMDFLAGS    = @SIMD_CFLAGS@
CPPFLAGS     = 
LDFLAGS      = 
DEFS         = -DHAVE_CONFIG_H
LIBS         = -lhmmer -leasel -ldivsufsort    -lm
MPILIBS      =  

AR           = /usr/bin/ar 
RANLIB       = ranlib
INSTALL      = /usr/bin/install -c

ESLDIR       = lib/easel
SADIR        = lib/libdivsufsort

MYLIBDIRS    = -L${top_builddir}/${ESLDIR} -L${top_builddir}/${SADIR} -L${top_builddir}/src
MYINCDIRS    = -I${top_builddir}/${ESLDIR} \
               -I${top_builddir}/${SADIR} \
	       -I${top_builddir}/src \
               -I${top_srcdir}/${ESLDIR} \
	       -I${top_srcdir}/${SADIR} \
	       -I${top_srcdir}/src 
MYLIBDEPS    = ${top_builddir}/src/libhmmer.a \
               ${top_builddir}/${ESLDIR}/libeasel.a \
               ${top_builddir}/${SADIR}/libdivsufsort.a

QUIET_SUBDIR0 = +${MAKE} -C #space separator after -c
QUIET_SUBDIR1 = 
ifndef V
	QUIET_CC      = @echo '    ' CC $@;
	QUIET_GEN     = @echo '    ' GEN $@;
	QUIET_AR      = @echo '    ' AR $@;
	QUIET_SUBDIR0 = +@subdir=
	QUIET_SUBDIR1 = ; echo '    ' SUBDIR $$subdir; \
		        ${MAKE} -C $$subdir
endif

.PHONY: all distclean clean tags-append

all:    ${PROGS}

.c.o: 
	${QUIET_CC}${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${DEFS} ${PTHREAD_CFLAGS} ${MYINCDIRS} -o $@ -c $<		

${PROGS}: %: %.o  ${MYLIBDEPS}
	${QUIET_GEN}${CC} ${CFLAGS} ${SIMDFLAGS} ${DEFS} ${LDFLAGS} ${MYLIBDIRS} -o $@ $@.o ${MPILIBS} ${LIBS}

distclean: clean
	-rm -f Makefile

clean:
	-rm -f *.o *~ ${PROGS} 
	-rm -f *.gcno

tags-append:
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/x-*
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/*.pl
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/*.c
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/*.in

################################################################
# @LICENSE@
# 
# SVN $URL$
# SVN $Id$
################################################################

