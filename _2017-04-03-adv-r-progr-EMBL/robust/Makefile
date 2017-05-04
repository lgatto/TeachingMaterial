%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

%.html: %.Rmd
	Rscript -e "rmarkdown::render('$^', output_format=rmarkdown::ioslides_presentation())"
	# Rscript -e "rmarkdown::render('$^', output_format=rmarkdown::html_document())"

%.R: %.Rmd
	Rscript -e "knitr::purl('$^')"
	perl -pi -e 's/^## (--.+)?//' $@

all: 
	make 01-debug.md 02-profile.md 03-funprog.md
	make 01-debug.R 02-profile.R 03-funprog.R 

.PHONY: all
