---
layout: page
title: Minimal make
tagline: A minimal tutorial on make
---

I would argue that the most important tool for reproducible research
is not [Sweave](http://www.stat.uni-muenchen.de/~leisch/Sweave/) or 
[knitr](http://yihui.name/knitr/) but
*[make](http://www.gnu.org/software/make)* (cf
[rake](http://rake.rubyforge.org)).

Consider, for example, all of the files associated with a manuscript.
In the simplest case, I would have an [R](http://r-project.org)
scripts for each figure plus a [LaTeX](http://www.latex-project.org)
file for the main text.  And then a [BibTeX](http://www.bibtex.org)
file.

Compiling the final PDF is a bit of work: 

- Run each R script through R to produce the relevant
- Run latex and then bibtex and then latex a couple of more times.

And the R scripts need to be run before latex is, and only if they've
changed.

### A simple example

[make](http://www.gnu.org/software/make) makes this easy.  In your
directory for the manuscript, you create a text file called `Makefile`
that looks something like the following (here using
[pdflatex](http://www.tug.org/applications/pdftex/)).

    mypaper.pdf: mypaper.bib mypaper.tex Figs/fig1.pdf Figs/fig2.pdf
        latex mypaper
        bibtex mypaper
        latex mypaper
        latex mypaper

    Figs/fig1.pdf: R/fig1.R
    	cd R;R CMD BATCH fig1.R fig1.Rout

    Figs/fig2.pdf: R/fig2.R
    	cd R;R CMD BATCH fig2.R fig2.Rout
