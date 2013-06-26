LATEXFILES = *.aux\
	*.bbl\
	*.blg\
	*.ilg\
	*.log\
	*.nlo\
	*.nls\
	*.toc\
	*.aux\
	R-programming.out\
	Rplots.pdf\
	*.dvi\
	*.map\
	*.figlist\
	*.dep\
	*.dpth\
	*.makefile\
	R-programming.tex

all: 
	make vig
	make r

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif

vig: R-programming.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('R-programming.Rnw');"
	## bibtex R-programming
	## "$(R_HOME)/bin/R" CMD pdflatex R-programming.tex
	"$(R_HOME)/bin/R" CMD pdflatex R-programming.tex
r:
	"$(R_HOME)/bin/R" --vanilla -e "knitr::purl('R-programming.Rnw')"

.PHONY: clean allclean 

clean:	
	rm -f $(LATEXFILES)
	rm -f .Rhistory
	rm -rf figure
	rm -f *~
	rm -f */*~
