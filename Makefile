%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

%.html: %.Rmd
	Rscript -e "rmarkdown::render('$^', output_format=rmarkdown::ioslides_presentation())"
	# Rscript -e "rmarkdown::render('$^', output_format=rmarkdown::html_document())"

all: 
	make 01-intro.md 02-funprog.md 03-debug.md 04-perf.md  unittesting.md rc.md deferred-eval.Rmd

.PHONY: all
