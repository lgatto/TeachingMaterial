LATEXFILES = *.aux\
	*.bbl\
	*.blg\
	*.ilg\
	*.log\
	*.nlo\
	*.nls\
	*.toc\
	*.aux\
	Rplots.pdf\
	*.dvi\
	*.map\
	*.out\
	*.figlist\
	*.dep\
	*.dpth\
	*.makefile\
	S4-tutorial.tex

all: 
	make s4

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

s4: 
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('S4-tutorial.Rnw');"
	## bibtex S4-tutorial
	## "$(R_HOME)/bin/R" CMD pdflatex S4-tutorial.tex
	"$(R_HOME)/bin/R" CMD pdflatex S4-tutorial.tex


.PHONY: clean allclean 

clean:	
	rm -f $(LATEXFILES)
	rm -f .Rhistory
	rm -rf figure
	rm -f *~
	rm -f */*~
