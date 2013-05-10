all: estimatek estimate clean

setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif


estimatek:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('estimatek.Rnw');"
	pdflatex estimatek.tex

estimate:
	"$(R_HOME)/bin/R" CMD Sweave estimate.Rnw
	pdflatex estimate.tex

clean:
	rm -f *.dvi *.aux *.log *.nav *.out *.snm *.toc *~ *vrb 
	rm -f estimate.tex estimatek.tex
	rm -f estimate-0*pdf
	rm -rf figure

allclean: clean
	rm -f estimate.pdf estimatek.pdf

