all:
	make 01-rstats.html 02-rstats.html 03-rstats.html README.html

%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

%.html: %.md
	Rscript -e "rmarkdown::render('$^')"


.PHONY: all
