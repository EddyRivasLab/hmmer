# All Makefile.in's in src/ subdirs include the identical
# block of makefile code below...

SHELL        = /bin/sh

prefix       = /usr/local
exec_prefix  = ${prefix}
datarootdir  = ${prefix}/share
bindir       = ${exec_prefix}/bin
libdir       = ${exec_prefix}/lib
includedir   = ${prefix}/include

CC           = gcc

CFLAGS       = -O3 -pthread -fPIC
#add -mavx and -mavx2 to enable avx instructions
#add -mavx512bw and -mavx512dq for AVX-512 This requires GCC 5 because they cleverly changed their AVX 512
#compiler flags between GCC 4.9 and 5.

#code may crash if you enable instructions that your hardware doesn't have, even if you don't use any of those
#instructions.  Best guess is that this is because enabling the instructions causes the compiler to worry
#about saving/restoring the appropriate type of vector registers
SIMDFLAGS    = -msse2 -msse3 
#-mavx -mavx2
#-mavx512bw -mavx512dq
CPPFLAGS     = 
LDFLAGS      = 
DEFS         = -DHAVE_CONFIG_H
LIBS         = ${LIBTMP} -lhmmer -leasel -ldivsufsort   -lm

AR           = /usr/bin/ar 
RANLIB       = ranlib
INSTALL      = /usr/bin/install -c

# in MYINCDIRs, we have build dirs first, because headers generated by
# ./configure (p7_config.h, esl_config.h) are in build dirs, not
# source dirs; and we want to find these first (in case an errant one is
# in the source tree). All other headers are in the source tree.
# 
ESLDIR       = lib/easel
SADIR        = lib/libdivsufsort
MYLIBDIRS    = -L${top_builddir}/${ESLDIR}\
               -L${top_builddir}/${SADIR}\
               ${LIBTMPDIR} \
               -L${top_builddir}/src
MYINCDIRS    = -I${top_builddir}/${ESLDIR} \
               -I${top_builddir}/${SADIR} \
	       -I${top_builddir}/src \
	       -I${top_srcdir}/${ESLDIR} \
	       -I${top_srcdir}/${SADIR} \
	       -I${top_srcdir}/src
MYLIBDEPS    = ${LIBTMPDEP} ${top_builddir}/src/libhmmer.a ${top_builddir}/${ESLDIR}/libeasel.a ${top_builddir}/${SADIR}/libdivsufsort.a


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


.PHONY: all dev check tests install uninstall distclean clean tags-append
.FORCE:

all:   libhmmer-${MODULE}.stamp ${LIBTMPDEP} ${PROGS}
dev:   libhmmer-${MODULE}.stamp ${LIBTMPDEP} ${PROGS} ${UTESTS} ${STATS} ${BENCHMARKS} ${EXAMPLES}
check: libhmmer-${MODULE}.stamp ${LIBTMPDEP} ${PROGS} ${UTESTS}
tests: libhmmer-${MODULE}.stamp ${LIBTMPDEP} ${PROGS} ${UTESTS}

# In some subdirs, ${LIBOBJS} is empty. We can't give an empty arg string to 'ar'.
# Hence the ! -z shell test on ${LIBOBJS} below.
libhmmer-${MODULE}.stamp: ${LIBOBJS}
	${QUIET}+if [ ! -z "${LIBOBJS}" ]; then \
	   echo '    ' AR $@ ;\
	   ${AR} -r ${top_builddir}/src/libhmmer.a $? > /dev/null 2>&1 ;\
	   ${RANLIB} ${top_builddir}/src/libhmmer.a ;\
	   echo ${MODULE} "objects compiled for libhmmer:\c" > $@ ;\
	   date >> $@ ;\
	else \
	   echo ${MODULE} "has no libhmmer objects\c" > $@ ;\
	fi


# The way we build unit tests, benchmarks, etc., our compile lines end
# up with 'cc -o foo_utest foo.c ${OBJS}' and ${OBJS} includes foo.o;
# compiler may barf on duplicate symbols in foo.c, foo.o. So, let the
# linker deal with the issue; linker should only seek for symbols
# that haven't already been defined. 
# When ${OBJS} is non-empty, the including Makefile also sets:
#    ${LIBTMP}    = -l./lib${MODULE}.a     prepended to ${LIBS}
#    ${LIBTMPDIR} = -L.                    added to ${MYLIBDIRS}
#    ${LIBTMPDEP} = lib${MODULE}.a         prepended to ${MYLIBDEPS}
lib${MODULE}.a: ${OBJS}
	${QUIET_AR}${AR} -r lib${MODULE}.a $? > /dev/null 2>&1 
	@${RANLIB} lib${MODULE}.a 

${LIBOBJS}:  ${HDRS}
${OBJS}:     ${HDRS}

.c.o:  
	${QUIET_CC}${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${DEFS} ${PTHREAD_CFLAGS} ${MYINCDIRS} -o $@ -c $<

${PROGS}: %: %.o  libhmmer-${MODULE}.stamp ${MYLIBDEPS}
	${QUIET_GEN}${CC} ${CFLAGS} ${SIMDFLAGS} ${DEFS} ${LDFLAGS} ${MYLIBDIRS} -o $@ $@.o ${LIBS}

${UTESTS}: libhmmer-${MODULE}.stamp ${MYLIBDEPS}
	@BASENAME=`echo $@ | sed -e 's/_utest//'| sed -e 's/^p7_//'` ;\
	DFLAG=`echo $${BASENAME} | sed -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`;\
	DFLAG=p7$${DFLAG}_TESTDRIVE ;\
	if test -e ${srcdir}/p7_$${BASENAME}.c; then \
           DFILE=${srcdir}/p7_$${BASENAME}.c ;\
        else \
           DFILE=${srcdir}/$${BASENAME}.c ;\
	fi;\
	if test ${V} ;\
	   then echo "${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS}" ;\
	   else echo '    ' GEN $@ ;\
	fi ;\
	${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS}

${STATS}: libhmmer-${MODULE}.stamp ${MYLIBDEPS}
	@BASENAME=`echo $@ | sed -e 's/_stats//' | sed -e 's/^p7_//'`;\
	DFLAG=`echo $${BASENAME} | sed -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`;\
	DFLAG=p7$${DFLAG}_STATS ;\
	if test -e ${srcdir}/p7_$${BASENAME}.c; then \
           DFILE=${srcdir}/p7_$${BASENAME}.c ;\
        else \
           DFILE=${srcdir}/$${BASENAME}.c ;\
	fi;\
	if test ${V} ;\
	   then echo "${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS}" ;\
	   else echo '    ' GEN $@ ;\
	fi ;\
	${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS}

${BENCHMARKS}: libhmmer-${MODULE}.stamp  ${MYLIBDEPS}
	@BASENAME=`echo $@ | sed -e 's/_benchmark//' | sed -e 's/^p7_//'`;\
	DFLAG=`echo $${BASENAME} | sed -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`;\
	DFLAG=p7$${DFLAG}_BENCHMARK ;\
	if test -e ${srcdir}/p7_$${BASENAME}.c; then \
           DFILE=${srcdir}/p7_$${BASENAME}.c ;\
        else \
           DFILE=${srcdir}/$${BASENAME}.c ;\
	fi;\
	if test ${V} ;\
	   then echo "${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS}" ;\
	   else echo '    ' GEN $@ ;\
	fi ;\
	${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS}

${EXAMPLES}: libhmmer-${MODULE}.stamp ${MYLIBDEPS}
	@BASENAME=`echo $@ | sed -e 's/_example[0-9]*//'| sed -e 's/^p7_//'` ;\
	DFLAG=`echo $${BASENAME} | sed -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`;\
	DFLAG=p7$${DFLAG}_EXAMPLE ;\
	if test -e ${srcdir}/p7_$${BASENAME}.c; then \
           DFILE=${srcdir}/p7_$${BASENAME}.c ;\
        else \
           DFILE=${srcdir}/$${BASENAME}.c ;\
	fi;\
	if test ${V} ;\
	   then echo "${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS}" ;\
	   else echo '    ' GEN $@ ;\
	fi ;\
	${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS}

install:
	${QUIET}if [ ! -z "${PROGS}" ]; then \
	   for file in ${PROGS}; do \
	      echo '    ' INSTALL $$file ;\
	      ${INSTALL} -m 0755 $$file ${DESTDIR}${bindir}/ ;\
	   done ;\
	fi

uninstall:
	${QUIET}if [ ! -z "${PROGS}" ]; then \
	   for file in ${PROGS}; do \
	      echo '    ' UNINSTALL $$file ;\
	      rm -f ${DESTDIR}${bindir}/$$file ;\
	   done ;\
	fi

distclean: clean
	-rm -f Makefile 

clean:
	-rm -f libhmmer-${MODULE}.stamp lib${MODULE}.a ${PROGS} ${UTESTS} ${STATS} ${BENCHMARKS} ${EXAMPLES}
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out cscope.out *.gcno *.gcda *.gcov
	${QUIET}for prog in ${PROGS} ${UTESTS} ${STATS} ${BENCHMARKS} ${EXAMPLES}; do \
	   if test -d $$prog.dSYM; then rm -rf $$prog.dSYM; fi ;\
	done

tags-append:
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/*.c ${srcdir}/*.h ${srcdir}/*.in


################################################################
# @LICENSE@
#
# SVN $URL$
# SVN $Id$
################################################################