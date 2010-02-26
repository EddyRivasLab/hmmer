################################################################
# Makefile for the HMMER testsuite
# Makefile is generated from Makefile.in by autoconf
# DO NOT EDIT Makefile: only edit Makefile.in
#
# SRE, Sat Feb 16 10:18:40 2008 [janelia]
# SVN $Id: Makefile.in 3152 2010-02-07 22:55:22Z eddys $
#

# VPATH and shell configuration
top_srcdir     = ..
srcdir         = .
VPATH          =  ..

SRCDIR   = ../src
ESLDIR   = ../easel
IMPLDIR  = ${SRCDIR}/impl_sse
SQCLEVEL = 2

PROGS    = 

ALIS     = \
	20aa\
	Caudal_act\
	LuxC\
	M1\
	Patched\
	RRM_1\
	SMC_N\
	XYPPX\
	ecori

all:  
	cd ${ESLDIR};  ${MAKE} easel_lib
	cd ${ESLDIR};  ${MAKE} easel_progs
	cd ${IMPLDIR}; ${MAKE} impl_lib
	cd ${SRCDIR};  ${MAKE} hmmer_lib
	cd ${SRCDIR};  ${MAKE} hmmer_progs
	cd ${SRCDIR};  ${MAKE} hmmer_utests
	cd ${IMPLDIR}; ${MAKE} impl_utests
	${MAKE} testsuite_progs

check: 
	cd ${ESLDIR};  ${MAKE} easel_lib
	cd ${ESLDIR};  ${MAKE} easel_progs
	cd ${ESLDIR};  ${MAKE} easel_utests
	cd ${IMPLDIR}; ${MAKE} impl_lib
	cd ${SRCDIR};  ${MAKE} hmmer_lib
	cd ${SRCDIR};  ${MAKE} hmmer_progs
	cd ${SRCDIR};  ${MAKE} hmmer_utests
	cd ${IMPLDIR}; ${MAKE} impl_utests
	${MAKE} testsuite_progs
	${MAKE} run_sqc

testsuite_progs: ${PROGS}

run_sqc:
	@echo 
	@echo Running HMMER test suite...
	@echo
	@${srcdir}/../easel/devkit/sqc ${SQCLEVEL} ${srcdir}/testsuite.sqc .. ${srcdir}/..

hmms:
	for ali in ${ALIS}; do\
	   ${top_srcdir}/src/hmmbuild $$ali.hmm $$ali.sto;\
	done

clean:
	-rm -f *.gcno *.gcda *.gcov
	-rm -f *.o *~ 

distclean:
	make clean
	-rm -f Makefile
