R_OPTS=--vanilla

mypaper.pdf: mypaper.bib mypaper.tex Figs/fig1.pdf Figs/fig2.pdf
	pdflatex mypaper
	bibtex mypaper
	pdflatex mypaper
	pdflatex mypaper

Figs/%.pdf: R/%.R
	cd $(<D);R CMD BATCH ${R_OPTS} $(<F)


