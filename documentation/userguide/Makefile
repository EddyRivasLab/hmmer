# Makefile for producing the HMMER documentation
# RCS $Id$
##########
# HMMER - Biological sequence analysis with profile HMMs
# Copyright (C) 1992-1998 Washington University School of Medicine
#
#   This source code is distributed under the terms of the 
#   GNU General Public License. See the files COPYING and 
#   LICENSE for details.
#    
###########

# Preparation notes.
#   - Install latex2html (version 98-1p1)
#   - latex2html.config is customized for HMMER User Guide (beware)
#   - Install PolyglotMan (version 3.0.5)
#   - need a copy of wulogo.gif for putting in HTML directory.



SHELL       = /bin/sh

MANS = hmmer hmmalign hmmbuild hmmcalibrate hmmconvert hmmemit hmmpfam hmmsearch 

all:
	make manpages
	make dvi
	make postscript
	make html

# manpages: convert man pages to LaTeX chapter in User Guide.
#    uses PolyglotMan 3.0.5 "rman", and rmanprocess.pl script
manpages:
				# write the manpages.tex "chapter" by conversion
	@echo "%% This chapter automatically generated. Do not edit." > manpages.tex
	@echo "\chapter{Manual pages}" >> manpages.tex 
	for file in $(MANS); do\
	   rman	-f LaTeX2e ../Man/$$file.man | ./rmanprocess.pl >> manpages.tex;\
	done

# dvi: make the dvi file
#
dvi:
	latex main
	bibtex main
	latex main
	latex main

# postscript: make the PostScript copy, Userguide.ps
#             "make manpages" and "make dvi" first
#
postscript:
	dvips -o Userguide.ps main


# html: invoke the LaTeX2html converter to create the hmmer-html/ HTML directory.
#       Do "make manpages" and "make dvi" first.
#
html:
	latex2html -split 4 -local_icons -t "HMMER 2.0 User's Guide" -dir hmmer-html main
	-cp wulogo.gif hmmer-html/
	tar cf hmmer-html.tar hmmer-html
	compress hmmer-html.tar

# clean: clean up after yourself
#
clean:
	-rm -f *.dvi *.aux *.toc *.bbl *.blg *.log *~
	-rm Userguide.ps
	-rm hmmer-html.tar.Z
	-rm -rf hmmer-html