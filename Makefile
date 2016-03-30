%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

slides.html: slides.md
	Rscript -e "rmarkdown::render('slides.Rmd')"

