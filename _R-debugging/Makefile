vectorisation.pdf: debugging.Rnw
	R --vanilla -e "library(knitr); knit('debugging.Rnw')"
	pdflatex debugging.tex

unittesting.md: unittesting.Rmd
		R --vanilla -e "library(knitr); knit2html('unittesting.Rmd')"

.PHONY: clean

clean:
	rm -f *~
	rm -f debugging.aux debugging.log debugging.nav debugging.org debugging.out debugging.snm debugging.tex debugging.toc debugging.vrb

allclean:
	make clean
	rm -rf cache figure

