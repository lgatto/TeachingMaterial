code.md: code.Rmd
		R --vanilla -e "library(knitr); knit('code.Rmd')"

all:
	make code.md

clean:
	rm -fr cache
	rm -f rprof
	rm -f *~
