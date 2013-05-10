setvars:
ifeq (${R_HOME},)
R_HOME= $(shell R RHOME)
endif


all: spr spr-4up.pdf fib fibonacci-4up.pdf spr.R fib.R rr

spr:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('spr.Rnw');"
	pdflatex spr.tex

fib:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('fibonacci.Rnw');"
	pdflatex fibonacci.tex

%-4up.pdf: %.pdf
	pdfnup -q --nup 2x2 --suffix '4up' $<

spr.R:
	R --vanilla -e 'Stangle("spr.Rnw")'

fib.R:
	R --vanilla -e 'Stangle("fibonacci.Rnw")'

clean:
	rm -f *.dvi *.aux *.log *.nav *.out *.snm *.toc *~ *vrb 
	rm -f */*~
	rm -rf figure cache
	rm -rf hist.pdf hist.png
	rm -rf spr.tex rr.tex fibonacci.tex

