%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

%.html: %.Rmd
	Rscript -e "rmarkdown::render('$^', output_format=rmarkdown::ioslides_presentation())"
	# Rscript -e "rmarkdown::render('$^', output_format=rmarkdown::html_document())"

%.R: %.Rmd
	Rscript -e "knitr::purl('$^')"
	perl -pi -e 's/^## (--.+)?//' $@

all: 
	make 01-intro.md 02-funprog.md 03-debug.md 04-perf.md
	make 01-intro.R 02-funprog.R 03-debug.R 04-perf.R 

.PHONY: all
