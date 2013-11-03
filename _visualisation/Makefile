setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif


ggplot2: ggplot2.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('ggplot2.Rnw');"
	pdflatex ggplot2.tex

clean:
	rm -f *.dvi *.aux *.log *.nav *.out *.snm *.toc *~ *vrb 
	rm -f */*~
	rm -rf figure cache


