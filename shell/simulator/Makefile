report.html: report.Rmd sim1.png sim2.png
	Rscript -e "require(knitr); knit2html('report.Rmd')"

sim1.dat: params.R simulator.R
	Rscript simulator.R rnorm > sim1.dat

sim2.dat: params.R simulator.R
	Rscript simulator.R runif > sim2.dat

sim1.png: sim1.dat plotter.R
	Rscript plotter.R sim1.dat

sim2.png: sim2.dat plotter.R
	Rscript plotter.R sim2.dat


## Old version for making PDF; needs targets for 
## making sim1.pdf from sim1.dat.
##report.pdf: report.tex sim1.pdf sim2.pdf 
##	texi2pdf report.tex



.PHONY: all clean

all: report.html

clean:
	rm -f report.html report.md
	rm -f report.pdf report.log report.aux
	rm -f sim1.* sim2*
