all:
	make rccpp.pdf
	make rc.md
	make deferred-eval.md

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

%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

rccpp.pdf: rccpp.Rnw intro.tex call.tex rcpp.tex
	R --vanilla -e "library(knitr); knit2pdf('rccpp.Rnw');"
	pdflatex rccpp.tex

intro.tex: intro.Rnw
	R --vanilla -e "library(knitr); knit('intro.Rnw');"

call.tex: call.Rnw
	R --vanilla -e "library(knitr); knit('call.Rnw');"

rcpp.tex: rcpp.Rnw
	R --vanilla -e "library(knitr); knit('rcpp.Rnw');"


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
