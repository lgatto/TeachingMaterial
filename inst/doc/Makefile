LATEXFILES = *.aux\
        *.bbl\
        *.blg\
        *.ilg\
        *.log\
        *.nlo\
        *.nls\
        *.toc\
        *.aux\
	*.map\
	*.sh\
	sequences.tex\
        *.out


all: 
	##R CMD pgfsweave sequences.Rnw
	"$(R_HOME)/bin/R" CMD Sweave sequences.Rnw
	pdflatex sequences.tex
	pdflatex sequences.tex

setvars:
ifeq (${R_HOME},)
R_HOME=	$(shell R RHOME)
endif

clean:
	rm -f $(LATEXFILES)
	rm sequences-gccount.pdf
	rm -f *~

