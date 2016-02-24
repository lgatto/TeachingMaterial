bmark-prof-optim.md: bmark-prof-optim.Rmd
		R --vanilla -e "library(knitr); knit('bmark-prof-optim.Rmd')"

all:
	make bmark-prof-optim.md

clean:
	rm -fr cache
	rm -f rprof
	rm -f *~
