all:
	make book
	make code

book:
	/opt/Rpatched/lib64/R/bin/Rscript -e "bookdown::render_book('index.Rmd')"

code:
	/opt/Rpatched/lib64/R/bin/Rscript -e "sapply(list.files(pattern = 'Rmd'), knitr::purl)"

