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
	make RIntro
	## make StangleAll

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif


RIntro.pdf: RIntro.Rnw Sec-AboutR.tex Sec-IntroR.tex Sec-DataTypes1.tex Sec-DataTypes2.tex Sec-objects.tex Sec-Programming.tex Sec-Packages.tex Sec-UseCase.tex Sec-BiocCase.tex Sec-RBioc.tex
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('RIntro.Rnw');" 
	pdflatex RIntro.tex ## to get rid of ~

Sec-AboutR.tex: Sec-AboutR.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-AboutR.Rnw');"

Sec-IntroR.tex: Sec-IntroR.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-IntroR.Rnw');"

Sec-DataTypes1.tex: Sec-DataTypes1.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-DataTypes1.Rnw');"

Sec-DataTypes2.tex: Sec-DataTypes2.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-DataTypes2.Rnw');"

Sec-objects.tex: Sec-objects.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-objects.Rnw');"

Sec-Programming.tex: Sec-Programming.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-Programming.Rnw');"

Sec-RBioc.tex: Sec-RBioc.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-RBioc.Rnw');"

Sec-Packages.tex: Sec-Packages.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-Packages.Rnw');"

Sec-UseCase.tex: Sec-UseCase.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-UseCase.Rnw');"

Sec-BiocCase.tex: Sec-BiocCase.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-BiocCase.Rnw');"

Sec-Plotting.tex: Sec-Plotting.Rnw
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit('Sec-Plotting.Rnw');"

clean:
	rm -f $(LATEXFILES)
	rm -f *~
	rm -rf figure
	rm -rf Exercises/aqmReport
	rm -rf Exercises/*~
	rm -rf Exercises/*pdf
	rm -rf Exercises/*rda
	rm -rf Exercises/figure


# StangleAll:
# 	"$(R_HOME)/bin/R" --vanilla -e "sapply(dir(pattern = 'Rnw'), Stangle)"
# 	cat Sec-IntroR.R Sec-DataTypes1.R Sec-DataTypes2.R Sec-objects.R Sec-UseCase.R Sec-Plotting.R Sec-Programming.R Sec-Packages.R Sec-RBioc.R > R-Basics_LG_2012.R
