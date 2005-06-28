# Makefile for producing the HMMER documentation
# CVS $Id$
##########
# @LICENSE@   
###########

# Preparation notes.
#   - need latex2html (version 98-1p1 or later)
#   - need PolyglotMan (version 3.0.5 or later)

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

# Currently, only the PDF version is distributed w/ HMMER.
#
all:    pdf


# pdf: make the PDF for distribution
#
pdf:  manpages
	pdflatex main
	bibtex main
	pdflatex main
	pdflatex main
	mv main.pdf Userguide.pdf

# manpages: convert man pages to LaTeX chapter in User Guide.
#    uses PolyglotMan 3.0.5 "rman", and rmanprocess.pl script
manpages:
				# write the manpages.tex "chapter" by conversion
	@echo "%% This chapter automatically generated. Do not edit." > manpages.tex
	@echo "\section{Manual pages}" >> manpages.tex 
	for file in $(HMMERMANS); do\
	   rman	-f LaTeX2e ../man/$$file.man | ../../ssdk/rmanprocess.pl >> manpages.tex;\
	done


# dvi: make the dvi file, main.dvi
#
dvi:  manpages
	latex main
	bibtex main
	latex main
	latex main

# postscript: make the PostScript copy, Userguide.ps
#             "make manpages" and "make dvi" first
#
postscript:   manpages dvi
	dvips -o Userguide.ps main


# html: invoke the LaTeX2html converter to create the hmmer-html/ HTML directory.
#       Do "make manpages" and "make dvi" first.
#
# Yes, you do need the latex and bibtex commands to run first;
# that generates .aux, etc files that are needed for bibliography,
# crossrefs, etc.
#
# Known problems:
#   - not processing list of FAQ's
#   - in manpages, --foo options are becoming -foo
#   - generated images are coming out w/ grey backgrounds
#   - \pageref crossrefs are making boxes instead of page numbers
#   - font converts to italic after first FAQ, and never comes back
#   - man pages should be segmented into different html files
#   - [rman -f HTML does a much better job formatting man pages, as opposed
#      to the two step process of rman -f LaTeX; latex2html!]
#
html:   manpages dvi
	latex2html -white -no_math -html_version 3.2,math -debug -split 4 -local_icons -t "HMMER 2.3 User's Guide" -dir hmmer-html main
	-cp -f wulogo.gif hmmer-html/
	-rm -f hmmer-html.tar.Z
	tar cf hmmer-html.tar hmmer-html
	compress -f hmmer-html.tar
	-rm -f main-expanded.tex

# clean: clean up after yourself
#
clean:
	-rm -f *.dvi *.aux *.toc *.bbl *.blg *.log *.faq *~
	-rm manpages.tex
	-rm Userguide.ps Userguide.pdf
	-rm hmmer-html.tar.Z
	-rm -rf hmmer-html