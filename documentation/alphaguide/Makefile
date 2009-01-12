
all:    pdf

pdf:  
	pdflatex main
	bibtex main
	pdflatex main
	pdflatex main
	mv main.pdf Alphaguide.pdf

clean:
	-rm -f main.aux main.bbl main.blg main.log main.toc *~
	-rm Alphaguide.pdf
