%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

%.html: %.md
	Rscript -e "rmarkdown::render('$^')"

index.html: README.md
	make README.html
	mv README.html index.html

