%.md: %.Rmd
	Rscript -e 'require(knitr); knit("$^")'

# rprog.md: rprog.Rmd
# 	Rscript -e 'require(knitr); knit("rprog.Rmd")'

.PHONY: all clean

all: intro.md ressources.md rprog.md viz.md

clean:
	rm -f *~
