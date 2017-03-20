# Top level Makefile
# 
# On most systems, to build code you should only need:
#     % ./configure; make
# Optionally, you can run a test suite:
#     % make check
# And optionally, you can install everything more permanently:
#     % make install
#


# VPATH and shell configuration
#
top_srcdir     = .
srcdir         = .

SHELL          = /bin/sh

# location of easel
ESLDIR         = lib/easel

# location of libdivsufsort for suffix array creation
SADIR          = lib/libdivsufsort


# Package information
#
PACKAGE         = HMMER
PACKAGE_VERSION = 4.0dev
PACKAGE_TARNAME = hmmer

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

# Compiler configuration
#
CC        = mpicc
CFLAGS    = -O3 -g -pthread -fPIC
LDFLAGS   = -static 
SIMDFLAGS = @SIMD_CFLAGS@
CPPFLAGS  = 

# Other tools
#
AR        = /usr/bin/ar 
RANLIB    = ranlib
INSTALL   = /usr/bin/install -c

# beautification magic stolen from git 
#
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

.PHONY: all dev check pdf install uninstall clean distclean TAGS 

# all: Compile all documented executables.
#      (Excludes test programs.)
#
all: 
	${QUIET_SUBDIR0}${ESLDIR}     ${QUIET_SUBDIR1} all
	${QUIET_SUBDIR0}${SADIR}      ${QUIET_SUBDIR1} all
	${QUIET_SUBDIR0}src           ${QUIET_SUBDIR1} all
	${QUIET_SUBDIR0}benchmarks    ${QUIET_SUBDIR1} all

# dev: compile all executables, including drivers.
#
dev: 
	${QUIET_SUBDIR0}${ESLDIR}  ${QUIET_SUBDIR1} dev
	${QUIET_SUBDIR0}${SADIR}   ${QUIET_SUBDIR1} all
	${QUIET_SUBDIR0}src        ${QUIET_SUBDIR1} dev
	${QUIET_SUBDIR0}benchmarks ${QUIET_SUBDIR1} all

# tests: compile all test drivers for 'make check'
#
tests:
	${QUIET_SUBDIR0}${ESLDIR}  ${QUIET_SUBDIR1} tests
	${QUIET_SUBDIR0}src        ${QUIET_SUBDIR1} tests

# check: Run test suites.
#
check:
	${QUIET_SUBDIR0}${ESLDIR}  ${QUIET_SUBDIR1} tests
	${QUIET_SUBDIR0}${SADIR}   ${QUIET_SUBDIR1} all	
	${QUIET_SUBDIR0}src        ${QUIET_SUBDIR1} tests
	${QUIET_SUBDIR0}${ESLDIR}  ${QUIET_SUBDIR1} check
	${QUIET_SUBDIR0}testsuite  ${QUIET_SUBDIR1} check

# pdf: compile the User Guides.
#
pdf:
	${QUIET_SUBDIR0}documentation ${QUIET_SUBDIR1} pdf

# install: installs the binaries in ${bindir}/
#          When man pages are done, will install man pages in MANDIR/man1/  (e.g. if MANSUFFIX is 1)
#          Creates these directories if they don't exist.
#          Prefix those paths with ${DESTDIR} (rarely used, usually null;
#          may be set on a make command line when building contrib RPMs).
install: 
	${INSTALL} -d ${DESTDIR}${bindir}
	${INSTALL} -d ${DESTDIR}${man1dir}
	${INSTALL} -d ${DESTDIR}${pdfdir}
	${QUIET_SUBDIR0}src           ${QUIET_SUBDIR1} install
	${QUIET_SUBDIR0}documentation ${QUIET_SUBDIR1} install

# uninstall: Reverses the steps of "make install".
#
uninstall: 
	${QUIET_SUBDIR0}src           ${QUIET_SUBDIR1} uninstall
	${QUIET_SUBDIR0}documentation ${QUIET_SUBDIR1} uninstall

# "make clean" removes almost everything except configuration files.
#
clean:
	${QUIET_SUBDIR0}src           ${QUIET_SUBDIR1} clean
	${QUIET_SUBDIR0}benchmarks    ${QUIET_SUBDIR1} clean
	${QUIET_SUBDIR0}testsuite     ${QUIET_SUBDIR1} clean
	${QUIET_SUBDIR0}documentation ${QUIET_SUBDIR1} clean
	${QUIET_SUBDIR0}${ESLDIR}     ${QUIET_SUBDIR1} clean
	${QUIET_SUBDIR0}${SADIR}      ${QUIET_SUBDIR1} clean
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out

# "make distclean" leaves a pristine source distribution.
#
distclean:
	${QUIET_SUBDIR0}src           ${QUIET_SUBDIR1} distclean
	${QUIET_SUBDIR0}benchmarks    ${QUIET_SUBDIR1} distclean
	${QUIET_SUBDIR0}testsuite     ${QUIET_SUBDIR1} distclean
	${QUIET_SUBDIR0}documentation ${QUIET_SUBDIR1} distclean
	${QUIET_SUBDIR0}${ESLDIR}     ${QUIET_SUBDIR1} distclean
	${QUIET_SUBDIR0}${SADIR}      ${QUIET_SUBDIR1} distclean
	-rm config.log config.status
	-rm -rf autom4te.cache
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out
	-rm -f cscope.po.out cscope.out cscope.in.out cscope.files
	-rm Makefile
#Use 'ifneq' instead of 'test -e' because the '+@' in QUIET_SUBDIR0 can't
#be passed to the shell. Note that ifneq breaks if indented.
ifneq (,$(wildcard ./documentation/release-notes/LICENSE.sh))
	-rm -f documentation/release-notes/LICENSE.sh
endif

TAGS:
	./makeTAGS.sh


################################################################
# @LICENSE@
#
# SVN $URL$
# SVN $Id$
################################################################
