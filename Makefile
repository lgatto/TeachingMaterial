all:
	make 00-intro.html
	make 10-data.html
	make 20-uml.html
	make 30-sml.html
	make 99-more.html

%.md: %.Rmd
	/opt/Rpatched/lib64/R/bin/Rscript -e "knitr::knit('$^')"

%.html: %.md
	/opt/Rpatched/lib64/R/bin/Rscript -e "rmarkdown::render('$^')"

