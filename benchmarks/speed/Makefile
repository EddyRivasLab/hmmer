top_srcdir = ../..
srcdir     = .


.PHONY: all distclean clean tags-append

all:

distclean: clean
	-rm -f Makefile

clean:
	-rm -f *~ 

tags-append:
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/x-*
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/*.pl
	etags -o ${top_srcdir}/TAGS -a ${srcdir}/*.in
