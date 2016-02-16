%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

%.html: %.Rmd
	Rscript -e "rmarkdown::render('$^', output_format=rmarkdown::ioslides_presentation())"
	# Rscript -e "rmarkdown::render('$^', output_format=rmarkdown::html_document())"
