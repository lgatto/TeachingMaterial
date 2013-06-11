LATEXFILES = *.aux\
        *.bbl\
        *.vrb\
	*.nav\
	*.snm\
        *.blg\
        *.ilg\
        *.log\
        *.nlo\
        *.nls\
        *.toc\
        *.aux\
	.Rhistory\
	Rplots.pdf\
	*.tex\
	*.dvi\
	*.map\
        *.out\
	*.tikz\

all:
	make pdf
	make StangleAll

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif

pdf: R-Basics.tex
	pdflatex R-Basics.tex

R-Basics.tex: Sec-Intro.tex Sec-DataTypes.tex Sec-Programming.tex Sec-RBioc.tex Sec-Packages.tex Sec-Manip.tex Sec-Useful.tex Sec-Plotting.tex R-Basics.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('R-Basics.Rnw');"

Sec-Intro.tex:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-Intro.Rnw');"

Sec-DataTypes.tex:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-DataTypes.Rnw');"

Sec-Programming.tex:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-Programming.Rnw');"

Sec-RBioc.tex:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-RBioc.Rnw');"

Sec-Packages.tex:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-Packages.Rnw');"

Sec-Manip.tex:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-Manip.Rnw');"

Sec-Useful.tex:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-Useful.Rnw');"

Sec-Plotting.tex:
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-Plotting.Rnw');"

clean:
	rm -f $(LATEXFILES)
	rm -f *~
	rm -rf figure
	rm -f Data/data.rda

StangleAll:
	"$(R_HOME)/bin/R" --vanilla -e "sapply(dir(pattern = 'Rnw'), Stangle)"
	cat Sec-Intro.R Sec-DataTypes.R Sec-Manip.R Sec-Useful.R Sec-Plotting.R Sec-Programming.R Sec-Packages.R Sec-RBioc.R > R-Basics_LG_2012.R
