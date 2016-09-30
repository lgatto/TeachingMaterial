all:
	make advr.pdf

advr.tex: advr.Rnw
	make clean
	R-3.3 CMD Sweave advr.Rnw

advr.pdf: advr.tex
	pdflatex advr.tex
	pdflatex advr.tex


clean:
	rm -rf rprof
	rm -rf myRpackage
	rm -rf advr-myCode.pdf advr-plotmsg.pdf advr.tex
	rm -f *.aux
	rm -f *.map
	rm -f *.log
	rm -f *.nav
	rm -f *.out
	rm -f *.snm
	rm -f *.toc
	rm -f *.tex.backup
	rm -f *.dvi
	rm -f *.vrb
	rm -f *.bbl
	rm -f *.blg
	rm -f *~
	rm -f Rplots.pdf

.PHONY: clean all
