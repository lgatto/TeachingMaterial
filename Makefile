all:
	make rccpp.pdf
	make rc.md

LATEXFILES = *.aux\
	*.aux\
	*.map\
	*.log\
	*.nav\
	*.out\
	*.snm\
	*.toc\
	*.tex.backup\
	*.dvi\
	*.vrb\
	*.bbl\
	*.blg\
	*~\
	Rplots.pdf\
	*tex

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif



rccpp.pdf: rccpp.Rnw intro.tex call.tex rcpp.tex
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('rccpp.Rnw');"
	pdflatex rccpp.tex

intro.tex: intro.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('intro.Rnw');"

call.tex: call.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('call.Rnw');"

rcpp.tex: rcpp.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('rcpp.Rnw');"

rc.md: rc.Rmd
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('rc.Rmd')"


clean:
	rm -f $(LATEXFILES)
	rm -rf cache
	rm -rf figure
	rm -f .Rhistory
	rm -f src/*~
	rm -f src/*.o
	rm -f src/*.so
	rm -rf mypackage
	rm -f src/.Rhistory

.PHONY: clean 
