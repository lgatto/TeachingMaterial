functional-programming.pdf: functional-programming.Rnw
	R --vanilla -e "library(knitr); knit('functional-programming.Rnw')"
	pdflatex functional-programming.tex


.PHONY: clean

clean:
	rm -f *~
	rm -f functional-programming.aux functional-programming.log functional-programming.nav functional-programming.org functional-programming.out functional-programming.snm functional-programming.tex functional-programming.toc functional-programming.vrb

allclean:
	make clean
	rm -rf cache figure
