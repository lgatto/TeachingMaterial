
%.md: %.Rmd
	/opt/Rpatched/lib/R/bin/Rscript -e 'require("knitr"); knit("$^")'

all: 00-before-we-start.md  01-intro-to-R.md  02-starting-with-data.md  03-data-frames.md  04-dplyr.md  05-datastructures.md


clean:
	rm -f *~
	rm -rf .Rcache

.PHONY: clean all
