expdes-slides.pdf: 01-expdesign.Rmd
	Rscript -e "rmarkdown::render('01-expdesign.Rmd', output_file = 'expdes-slides.pdf', output_format = 'beamer_presentation')"

expdes-handouts.pdf: 01-expdesign.Rmd
	Rscript -e "rmarkdown::render('01-expdesign.Rmd', output_file = 'expdes-handouts.pdf', output_format = 'pdf_document')"

all: 
	make expdes-slides.pdf
	make expdes-handouts.pdf
