%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"
