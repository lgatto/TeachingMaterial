r-knitr.md: r-knitr.Rmd
	Rscript -e "knitr::knit('r-knitr.Rmd')"

r-knitr.pdf: r-knitr.md
	Rscript -e "rmarkdown::render('r-knitr.md', output_format = 'pdf_document')"

r-knitr.html: r-knitr.md
	Rscript -e "rmarkdown::render('r-knitr.md', output_format = 'html_document')"

clean:
	rm -f r-knitr.html r-knitr.pdf r-knitr.md Rplots.pdf
	rm -f 01-intro.pdf  02-open-science.pdf  03-rr.pdf  04-ccl.pdf  05-refs.pdf
	rm -rf figure
	rm -f *~

01-intro.pdf: 01-intro.md
	pandoc $^ -o $@

02-open-science.pdf: 02-open-science.md
	pandoc $^ -o $@

03-rr.pdf: 03-rr.md
	pandoc $^ -o $@

04-ccl.pdf: 04-ccl.md
	pandoc $^ -o $@

05-refs.pdf: 05-refs.md
	pandoc $^ -o $@

open-rr-best-practice.pdf: 01-intro.pdf 02-open-science.pdf 03-rr.pdf 04-ccl.pdf 05-refs.pdf
	pdftk 0*.pdf cat output open-rr-best-practice.pdf

all:
	make open-rr-best-practice.pdf
	make clean

.PHONY: clean all
