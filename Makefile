
%.md: %.Rmd
	/opt/Rpatched/lib/R/bin/Rscript -e 'require("knitr"); knit("$^")'

%.html: %.md
	/opt/Rpatched/lib/R/bin/Rscript -e 'rmarkdown::render("$^")'

rintro:
	cd R-intro; make all

r4p:
	cd r4p; make all

clean:
	rm -f *~
	rm -f */*~
	rm -rf .Rcache
	rm -f .Rhistory
	rm -f */.Rhistory

.PHONY: clean rintro
