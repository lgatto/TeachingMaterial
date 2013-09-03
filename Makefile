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
	rm -rf QuickPackage
	rm -rf QuickPackage.Rcheck
	make QuickPackage
	make QuickPackage-4up.pdf
	make QuickPackage.R
	make clean

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif

QuickPackage:
	rm -rf QuickPackage
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('QuickPackage.Rnw');"
	pdflatex QuickPackage.tex

QuickPackageAndMore:
	rm -rf QuickPackage
	"$(R_HOME)/bin/R" --vanilla -e "library(knitr); knit2pdf('QuickPackageAndMore.Rnw');"
	pdflatex QuickPackageAndMore.tex


QuickPackage.R:
	R --vanilla -e 'Stangle("QuickPackage.Rnw")'

%-4up.pdf: %.pdf
	pdfnup -q --nup 2x2 --suffix '4up' $<

clean:
	rm -rf QuickPackage.Rcheck
	rm -f $(LATEXFILES)
	rm -f *~
	rm -rf sweave-cache

allclean:
	make clean
	rm -f QuickPackage-4up.pdf  
	## rm -f QuickPackage.pdf
	rm -rf QuickPackage
