%.md: %.Rmd
	Rscript -e 'require(knitr); options(knitr.unnamed.chunk.label = "$^"); knit("$^")'

%.html: %.md
	Rscript -e 'require(knitr); options(knitr.unnamed.chunk.label = "$^"); knit2html("$^")'

.PHONY: all clean

all: intro.md rprog.md viz.md

clean:
	rm -f *~
	rm -f *txt
	rm -f *html
