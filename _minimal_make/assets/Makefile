mypaper.pdf: mypaper.bib mypaper.tex Figs/fig1.pdf Figs/fig2.pdf
	pdflatex mypaper
	bibtex mypaper
	pdflatex mypaper
	pdflatex mypaper

Figs/fig1.pdf: R/fig1.R
	cd R;R CMD BATCH fig1.R fig1.Rout

Figs/fig2.pdf: R/fig2.R
	cd R;R CMD BATCH fig2.R fig2.Rout


