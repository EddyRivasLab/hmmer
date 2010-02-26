# documentation/userguide/Makefile.  Generated from Makefile.in by configure.
#
# Makefile for user documentation
#
# SRE, Tue Mar 10 08:56:51 2009 [Janelia]
# SVN $Id$

# VPATH and shell configuration
top_srcdir     = ../..
srcdir         = .


HMMER_VERSION   =  3.0rc1
HMMER_DATE      = "February 2010"
HMMER_COPYRIGHT = "Copyright (C) 2010 Howard Hughes Medical Institute."

SEDITION      = ../../easel/devkit/sedition

HMMERMANS     =\
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

all:    pdf

pdf:   titlepage.tex copyright.tex manpages
	pdflatex main
	bibtex main
	pdflatex main
	pdflatex main
	mv main.pdf Userguide.pdf

titlepage.tex: titlepage.tex.in
	cp -f titlepage.tex.in titlepage.tex
	${SEDITION} HMMER_VERSION ${HMMER_VERSION} HMMER_DATE ${HMMER_DATE} titlepage.tex

copyright.tex: copyright.tex.in
	cp -f copyright.tex.in copyright.tex
	${SEDITION} HMMER_COPYRIGHT ${HMMER_COPYRIGHT} copyright.tex

# manpages: convert man pages to LaTeX chapter in User Guide.
#    uses PolyglotMan 3.2 "rman", and rmanprocess.pl script in easel's devkit
manpages:
	@echo "%% This chapter automatically generated. Do not edit." > manpages.tex
	@echo "\section{Manual pages}" >> manpages.tex 
	for file in $(HMMERMANS); do\
	   rman	-f LaTeX2e ../man/$$file.man | ../../easel/devkit/rmanprocess.pl >> manpages.tex;\
	done

clean:
	-rm -f main.aux main.bbl main.blg main.log main.toc *~
	-rm Userguide.pdf

distclean:
	-rm -f main.aux main.bbl main.blg main.log main.toc *~
	-rm Userguide.pdf
	-rm -f titlepage.tex copyright.tex manpages.tex
	-rm -f Makefile

TAGS:
	etags *.tex Makefile.in
