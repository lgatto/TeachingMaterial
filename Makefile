vectorisation.pdf: vectorisation.Rnw
	R --vanilla -e "library(knitr); knit('vectorisation.Rnw')"
	pdflatex vectorisation.tex


.PHONY: clean

clean:
	rm -f *~
	rm -f vectorisation.aux vectorisation.log vectorisation.nav vectorisation.org vectorisation.out vectorisation.snm vectorisation.tex vectorisation.toc vectorisation.vrb

allclean:
	make clean
	rm -rf cache figure
