# Makefile for HMMER src directory
#
# SRE, Mon Jan  1 19:05:39 2007 [Casa de Gatos]
# SVN $Id: Makefile.in 3266 2010-05-06 14:50:24Z castellanos $

top_srcdir = ..
srcdir     = .
VPATH      =  ..

CC       = gcc -std=gnu99
CFLAGS   = -g -O3 -fomit-frame-pointer -fstrict-aliasing -mpentiumpro 
SIMDFLAGS= -msse2
CPPFLAGS = 
LDFLAGS  = 
DEFS     = -DHAVE_CONFIG_H
LIBS     =    -lm

AR       = /usr/bin/ar rcv 
RANLIB   = ranlib

ESLDIR    = ../easel
IMPLDIR   = impl_sse

# Configuration for optional pthreads multiprocessor support
#
MPILIBS   = 

SHELL     = /bin/sh
MYLIBDIRS = -L${ESLDIR} -L./${IMPLDIR} -L.
MYINCDIRS = -I${ESLDIR} \
	    -I${top_srcdir}/easel \
	    -I. \
	    -I${srcdir}

LIBS  += -lhmmer -lhmmerimpl -leasel
LIBDEP = ./libhmmer.a impl_sse/libhmmerimpl.a ${ESLDIR}/libeasel.a

PROGS = hmmalign\
	hmmbuild\
	hmmconvert\
	hmmemit\
	hmmfetch\
	hmmpress\
	hmmscan\
	hmmsearch\
	hmmsim\
	hmmstat\
	jackhmmer\
	phmmer\
	shmmer

HDRS =  hmmer.h\
	p7_config.h

OBJS =  build.o\
	emit.o\
	errors.o\
	evalues.o\
	eweight.o\
	generic_decoding.o\
	generic_fwdback.o\
	generic_null2.o\
	generic_msv.o\
	generic_optacc.o\
	generic_stotrace.o\
	generic_viterbi.o\
	generic_vtrace.o\
	generic_smithwaterman.o\
	generic_miyazawa.o\
	h2_io.o\
	heatmap.o\
	hmmer.o\
	logsum.o\
	modelconfig.o\
	modelstats.o\
	mpisupport.o\
	rankorder.o\
	seqmodel.o\
	tracealign.o\
	p7_alidisplay.o\
	p7_bg.o\
	p7_builder.o\
	p7_domaindef.o\
	p7_gmx.o\
	p7_hmm.o\
	p7_hmmfile.o\
	p7_pipeline.o\
	p7_prior.o\
	p7_profile.o\
	p7_spensemble.o\
	p7_tophits.o\
	p7_trace.o
#	island.o\


STATS = \
	evalues_stats

BENCHMARKS = \
	evalues_benchmark\
	generic_decoding_benchmark\
	generic_fwdback_benchmark\
	generic_msv_benchmark\
	generic_null2_benchmark\
	generic_optacc_benchmark\
	generic_stotrace_benchmark\
	generic_viterbi_benchmark

UTESTS =\
	build_utest\
	generic_fwdback_utest\
	generic_miyazawa_utest\
	generic_msv_utest\
	generic_stotrace_utest\
	generic_viterbi_utest\
	hmmer_utest\
	logsum_utest\
	modelconfig_utest\
	p7_gmx_utest\
	p7_hmm_utest\
	p7_hmmfile_utest\
	p7_profile_utest\
	p7_tophits_utest\
	p7_trace_utest

ITESTS = \
	itest_brute

EXAMPLES = \
	build_example\
	generic_fwdback_example\
	generic_msv_example\
	generic_optacc_example\
	generic_stotrace_example\
	generic_viterbi_example\
	generic_vtrace_example


all: libraries programs

# check: Run test suites.
check:
			${MAKE} all
			${MAKE} hmmer_utests
			${MAKE} hmmer_itests
			${MAKE} testsuite
	cd ${TESTDIR};  ${MAKE} run_sqc

dev: libraries programs utests itests misc

libraries: easel_lib impl_lib hmmer_lib

easel_lib:
	cd ${ESLDIR}; ${MAKE} libeasel.a

impl_lib:
	cd ${IMPLDIR}; ${MAKE} impl_lib

hmmer_lib: libhmmer.a

programs: libraries ${PROGS}

utests: libraries impl_utests ${UTESTS}

itests: libraries             ${ITESTS}

misc: libraries ${STATS} impl_benchmarks ${BENCHMARKS} impl_examples ${EXAMPLES}

# the following rules are called from the upper level makefile and
# any library dependancies are assume to be resolved.
hmmer_progs:       ${PROGS}
hmmer_stats:       ${STATS}
hmmer_utests:      impl_utests     ${UTESTS}
hmmer_itests:                      ${ITESTS}
hmmer_benchmarks:  impl_benchmarks ${BENCHMARKS}
hmmer_examples:    impl_examples   ${EXAMPLES}

libhmmer.a: ${OBJS} 
	${AR} libhmmer.a $?
	${RANLIB} libhmmer.a
	chmod 644 libhmmer.a

impl_utests:
	cd ${IMPLDIR}; ${MAKE} impl_utests

impl_benchmarks:
	cd ${IMPLDIR}; ${MAKE} impl_benchmarks

impl_examples:
	cd ${IMPLDIR}; ${MAKE} impl_examples

${OBJS}: ${HDRS}

${PROGS}: %: %.o  ${LIBDEP}
	${CC} ${CFLAGS} ${SIMDFLAGS} ${DEFS} ${LDFLAGS} ${MYLIBDIRS} -o $@ $@.o ${MPILIBS} ${LIBS}

.c.o:
	${CC} ${CFLAGS} ${SIMDFLAGS} ${DEFS} ${CPPFLAGS} ${MYINCDIRS} -o $@ -c $<

.a:	libraries

${ITESTS}: %: %.o  ${LIBDEP}
	${CC} ${CFLAGS} ${SIMDFLAGS} ${DEFS} ${LDFLAGS} ${MYLIBDIRS} -o $@ $@.o ${MPILIBS} ${LIBS}

${UTESTS}: ${LIBDEP}
	@BASENAME=`echo $@ | sed -e 's/_utest//' | sed -e 's/^p7_//'` ;\
	DFLAG=`echo $${BASENAME} | sed -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`;\
	DFLAG=p7$${DFLAG}_TESTDRIVE ;\
	if test -e ${srcdir}/p7_$${BASENAME}.c; then \
           DFILE=${srcdir}/p7_$${BASENAME}.c ;\
        else \
           DFILE=${srcdir}/$${BASENAME}.c ;\
	fi;\
	echo ${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS};\
	${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS}

${BENCHMARKS}: ${LIBDEP}
	@BASENAME=`echo $@ | sed -e 's/_benchmark//' | sed -e 's/^p7_//'` ;\
	DFLAG=`echo $${BASENAME} | sed -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`;\
	DFLAG=p7$${DFLAG}_BENCHMARK ;\
	if test -e ${srcdir}/p7_$${BASENAME}.c; then \
           DFILE=${srcdir}/p7_$${BASENAME}.c ;\
        else \
           DFILE=${srcdir}/$${BASENAME}.c ;\
	fi;\
	echo ${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS};\
	${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS}

${STATS}: ${LIBDEP}
	@BASENAME=`echo $@ | sed -e 's/_stats//' | sed -e 's/^p7_//'` ;\
	DFLAG=`echo $${BASENAME} | sed -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`;\
	DFLAG=p7$${DFLAG}_STATS ;\
	if test -e ${srcdir}/p7_$${BASENAME}.c; then \
           DFILE=${srcdir}/p7_$${BASENAME}.c ;\
        else \
           DFILE=${srcdir}/$${BASENAME}.c ;\
	fi;\
	echo ${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS};\
	${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS}

${EXAMPLES}: ${LIBDEP}
	@BASENAME=`echo $@ | sed -e 's/_example//' | sed -e 's/^p7_//'` ;\
	DFLAG=`echo $${BASENAME} | sed -e 'y/abcdefghijklmnopqrstuvwxyz/ABCDEFGHIJKLMNOPQRSTUVWXYZ/'`;\
	DFLAG=p7$${DFLAG}_EXAMPLE ;\
	if test -e ${srcdir}/p7_$${BASENAME}.c; then \
           DFILE=${srcdir}/p7_$${BASENAME}.c ;\
        else \
           DFILE=${srcdir}/$${BASENAME}.c ;\
	fi;\
	echo ${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS};\
	${CC} ${CFLAGS} ${SIMDFLAGS} ${CPPFLAGS} ${LDFLAGS} ${DEFS} ${MYLIBDIRS} ${MYINCDIRS} -D$${DFLAG} -o $@ $${DFILE} ${LIBS}

distclean:
	-rm -f *.o *~ Makefile.bak core ${PROGS} TAGS gmon.out libhmmer.a
	-rm -f ${UTESTS}
	-rm -f ${ITESTS}
	-rm -f ${STATS}
	-rm -f ${BENCHMARKS}
	-rm -f ${EXAMPLES}
	-rm -f *.gcno
	-rm -f Makefile p7_config.h
	cd  ${IMPLDIR}; ${MAKE} distclean

clean:
	-rm -f *.o *~ Makefile.bak core ${PROGS} TAGS gmon.out libhmmer.a
	-rm -f ${UTESTS}
	-rm -f ${ITESTS}
	-rm -f ${STATS}
	-rm -f ${BENCHMARKS}
	-rm -f ${EXAMPLES}
	-rm -f *.gcno
	cd ${IMPLDIR}; ${MAKE} clean

binclean:
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out libhmmer.a
	cd ${IMPLDIR}; ${MAKE} binclean

TAGS:
	etags  *.c *.h Makefile.in
	etags -a impl_sse/*.c   impl_sse/*.h   impl_sse/Makefile.in
	etags -a impl_vmx/*.c   impl_vmx/*.h   impl_vmx/Makefile.in
	etags -a impl_dummy/*.c impl_dummy/*.h impl_dummy/Makefile.in

###########
# @LICENSE@
###########
