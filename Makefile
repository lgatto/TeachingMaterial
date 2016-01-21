expdes-slides.pdf: 01-expdesign.Rmd
	Rscript -e "rmarkdown::render('01-expdesign.Rmd', output_file = 'expdes-slides.pdf', output_format = 'beamer_presentation')"

expdes-handouts.pdf: 01-expdesign.Rmd
	Rscript -e "rmarkdown::render('01-expdesign.Rmd', output_file = 'expdes-handouts.pdf', output_format = 'pdf_document')"

test-slides.pdf: 02-testing.Rmd
	Rscript -e "rmarkdown::render('02-testing.Rmd', output_file = 'test-slides.pdf', output_format = 'beamer_presentation')"

test-handouts.pdf: 02-testing.Rmd
	Rscript -e "rmarkdown::render('02-testing.Rmd', output_file = 'test-handouts.pdf', output_format = 'pdf_document')"

all: 
	make expdes-slides.pdf
	make expdes-handouts.pdf
	make test-slides.pdf
	make test-handouts.pdf
