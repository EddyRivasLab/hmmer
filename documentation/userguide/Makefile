# Makefile for producing the HMMER documentation
# RCS $Id$
##########
# HMMER - Biological sequence analysis with profile HMMs
# Copyright (C) 1992-2001 Washington University School of Medicine
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
#   - Install ps2pdf and Ghostscript to generate PDF copy of Userguide


SHELL       = /bin/sh

HMMERMANS = hmmer\
	hmmalign\
	hmmbuild\
	hmmcalibrate\
	hmmconvert\
	hmmemit\
	hmmfetch\
	hmmindex\
	hmmpfam\
	hmmsearch

SQUIDMANS = afetch\
	alistat\
	seqstat\
	sfetch\
	shuffle\
	sindex\
	sreformat	

all:
	make manpages
	make dvi
	make postscript
	make html
	make pdf

# manpages: convert man pages to LaTeX chapter in User Guide.
#    uses PolyglotMan 3.0.5 "rman", and rmanprocess.pl script
manpages:
				# write the manpages.tex "chapter" by conversion
	@echo "%% This chapter automatically generated. Do not edit." > manpages.tex
	@echo "\chapter{Manual pages}" >> manpages.tex 
	for file in $(HMMERMANS); do\
	   rman	-f LaTeX2e ../man/$$file.man | ./rmanprocess.pl >> manpages.tex;\
	done
	for file in $(SQUIDMANS); do\
	   rman -f LaTeX2e ../../squid/Man/$$file.man | ./rmanprocess.pl >> manpages.tex;\
	done

# dvi: make the dvi file
#
dvi:
	latex main
	bibtex main
	latex main
	latex main

# pdf: make the PDF for distribution
#
pdf:
	make manpages
	pdflatex main
	bibtex main
	pdflatex main
	pdflatex main
	mv main.pdf Userguide.pdf

# postscript: make the PostScript copy, Userguide.ps
#             "make manpages" and "make dvi" first
#
postscript: 
	dvips -o Userguide.ps main


# html: invoke the LaTeX2html converter to create the hmmer-html/ HTML directory.
#       Do "make manpages" and "make dvi" first.
#
# There appears to be a bug in latex2html with \input{} statements
# and multifile LaTeX documents. My workaround is to swallow the \input's
# into a single large temporary LaTeX file using the latex-swallow-inputs
# script.
# 
html:
	latex-swallow-inputs main > main-expanded.tex
	latex main-expanded
	bibtex main-expanded
	latex main-expanded
	latex main-expanded
	latex2html -debug -split 4 -local_icons -t "HMMER 2.1 User's Guide" -dir hmmer-html main-expanded
	-cp -f wulogo.gif hmmer-html/
	-rm -f hmmer-html.tar.Z
	tar cf hmmer-html.tar hmmer-html
	compress -f hmmer-html.tar
	-rm -f main-expanded.tex

# clean: clean up after yourself
#
clean:
	-rm -f *.dvi *.aux *.toc *.bbl *.blg *.log *~
	-rm Userguide.ps Userguide.pdf
	-rm hmmer-html.tar.Z
	-rm -rf hmmer-html