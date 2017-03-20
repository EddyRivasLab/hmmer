top_srcdir = ../..
srcdir     = .

SHELL      = /bin/sh

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

INSTALL     = /usr/bin/install -c

MANS =  hmmer\
	hmmalign\
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
	phmmer

ifndef V
	QUIET         = @
endif

.PHONY:  install uninstall clean distclean

install:
	${QUIET}for file in ${MANS}; do \
	   echo '    ' INSTALL $$file ;\
	   ${INSTALL} -m 0755 ${srcdir}/$$file.man ${DESTDIR}${man1dir}/$${file}${man1ext} ;\
	done

uninstall:
	${QUIET}for file in ${MANS}; do \
	   echo '    ' UNINSTALL $$file ;\
	   rm ${DESTDIR}${man1dir}/$${file}${man1ext} ;\
	done

clean:
	-rm -f *~ 

distclean:
	-rm -f *~ 
	-rm -f Makefile

################################################################
# @LICENSE@
#
# SVN $URL$
# SVN $Id$
################################################################
