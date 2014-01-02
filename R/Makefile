%.md: %.Rmd
	Rscript -e 'require(knitr); knit("$^")'

%.html: %.md
	Rscript -e 'require(knitr); knit2html("$^")'

.PHONY: all clean

all: intro.md ressources.md rprog.md viz.md

clean:
	rm -f *~
	rm -f *txt
