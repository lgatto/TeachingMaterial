all:
	make rpd.pdf

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


rpd.pdf: rpd.Rnw dist.tex rd.tex struct.tex testing.tex
	rm -rf myRpackage
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('rpd.Rnw');"
	pdflatex rpd.tex

struct.tex: struct.Rnw
	rm -rf myRpackage
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('struct.Rnw');"

dist.tex: dist.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('dist.Rnw');"

testing.tex: testing.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('testing.Rnw');"

rd.tex: rd.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('rd.Rnw');"

clean:
	rm -f $(LATEXFILES)
	rm -rf myRpackage
	rm -rf figure
