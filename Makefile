
%.md: %.Rmd
	/opt/Rpatched/lib/R/bin/Rscript -e 'require("knitr"); knit("$^")'

clean:
	rm -f *~
	rm -rf .Rcache

.PHONY: clean all
