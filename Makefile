parallel.pdf: parallel.Rnw
	R --vanilla -e "library(knitr); knit('parallel.Rnw')"
	pdflatex parallel.tex

clean:
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
	rm -f *.tex
	rm -f *~
	rm -f Rplots.pdf
	rm -f .Rhistory
	rm -f src/.Rhistory
