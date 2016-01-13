r-knitr.md: r-knitr.Rmd
	Rscript -e "knitr::knit('r-knitr.Rmd')"

r-knitr.pdf: r-knitr.md
	Rscript -e "rmarkdown::render('r-knitr.md', output_format = 'pdf_document')"

r-knitr.html: r-knitr.md
	Rscript -e "rmarkdown::render('r-knitr.md', output_format = 'html_document')"

clean:
	rm -f r-knitr.html r-knitr.pdf r-knitr.md Rplots.pdf
	rm -rf figure
	rm -f *~
