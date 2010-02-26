# Makefile for HMMER profmark benchmark generation tools
# The Makefile is generated from Makefile.in by autoconf
# DO NOT EDIT Makefile; only edit Makefile.in
#
# SRE, Mon Mar 31 09:12:28 2008 [Janelia]
# SVN $Id: Makefile.in 3152 2010-02-07 22:55:22Z eddys $

top_srcdir = ..
srcdir     = .
VPATH      =  ..

CC       = gcc -std=gnu99
CFLAGS   = -O3 -fomit-frame-pointer -malign-double -fstrict-aliasing
CPPFLAGS = 
LDFLAGS  = 
SIMDFLAGS= -msse2
DEFS     = -DHAVE_CONFIG_H
LIBS     =   -lm
MPILIBS  = 
IMPLDIR  = impl_sse

SRCDIR   = ../src
ESLDIR   = ../easel

ESLINC   = -I${ESLDIR} -I${top_srcdir}/easel
SRCINC   = -I${SRCDIR} -I${top_srcdir}/src

PROGS    = create-profmark\
	   rocplot

all: ${PROGS} 
	cd ${ESLDIR};            ${MAKE} easel_lib
	cd ${SRCDIR}/${IMPLDIR}; ${MAKE} impl_lib
	cd ${SRCDIR};            ${MAKE} hmmer_lib
	cd ${SRCDIR};            ${MAKE} hmmer_progs
	${MAKE} profmark_progs

profmark_progs: ${PROGS} 

.c.o: 
	${CC} ${ESLINC} ${SRCINC} ${CFLAGS} ${SIMDFLAGS} ${DEFS} -o $@ -c $<		

${PROGS}: %: %.o  ${SRCDIR}/libhmmer.a ${ESLDIR}/libeasel.a
	${CC} ${CFLAGS} ${DEFS} ${LDFLAGS} -L${ESLDIR} -L${SRCDIR} -o $@ $@.o -lhmmer -leasel ${LIBS}

distclean:
	make clean
	-rm -f Makefile

clean:
	-rm -f *.o *~  ${PROGS} TAGS
	-rm -f *.gcno

