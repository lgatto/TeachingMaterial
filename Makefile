all:
	make roo.pdf

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


roo.pdf: roo.Rnw intro.tex revision.tex oo-intro.tex roo.Rnw S3.tex S4.tex Ref.tex 
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('roo.Rnw');"
	pdflatex roo.tex

intro.tex: intro.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('intro.Rnw');"

revision.tex: revision.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('revision.Rnw');"

oo-intro.tex: oo-intro.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('oo-intro.Rnw');"

S3.tex: S3.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('S3.Rnw');"

S4.tex: S4.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('S4.Rnw');"

Ref.tex: Ref.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Ref.Rnw');"

clean:
	rm -f $(LATEXFILES)
	rm -rf figure
