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
	make RIntro.pdf
	## make StangleAll

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif


RIntro.pdf: RIntro.Rnw Sec-AboutR.tex Sec-IntroR.tex Sec-DataTypes1.tex Sec-DataTypes2.tex Sec-objects.tex Sec-Programming.tex Sec-Packages.tex Sec-UseCase.tex Sec-BiocCase.tex Sec-RBioc.tex
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('RIntro.Rnw');" 
	pdflatex RIntro.tex ## to get rid of ~

%.tex: %.Rnw
	Rscript -e 'require(knitr);  knit("$^")'


.PHONY: all clean handouts

handouts:
	make all
	cp -f RIntro.pdf handouts/.
	cp -f Exercises/*html handouts/ex/.
	cp -f Exercises/*CEL handouts/ex/.
	cp -f Exercises/*tsv handouts/ex/.
	cp -f Exercises/*csv handouts/ex/.
	echo "Course web page: http://lgatto.github.io/RIntro/" > handouts/README.txt
	cp -r RefCards handouts/.
	zip -r handouts handouts/*
	make clean

clean:
	rm -f $(LATEXFILES)
	rm -f *~
	rm -rf figure
	rm -rf Exercises/aqmReport
	rm -rf Exercises/*~
	rm -rf Exercises/*pdf
	rm -rf Exercises/*rda
	rm -rf Exercises/.Rhistory

# StangleAll:
# 	"$(R_HOME)/bin/R" --vanilla -e "sapply(dir(pattern = 'Rnw'), Stangle)"
# 	cat Sec-IntroR.R Sec-DataTypes1.R Sec-DataTypes2.R Sec-objects.R Sec-UseCase.R Sec-Plotting.R Sec-Programming.R Sec-Packages.R Sec-RBioc.R > R-Basics_LG_2012.R
